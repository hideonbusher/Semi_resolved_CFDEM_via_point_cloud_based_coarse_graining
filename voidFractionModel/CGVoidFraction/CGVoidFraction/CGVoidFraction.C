/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "CGVoidFraction.H"
#include "addToRunTimeSelectionTable.H"
#include <sstream>
#include <string>
#include "dictionary.H"
#include "IOobject.H"
#include "List.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(CGVoidFraction, 0);

addToRunTimeSelectionTable
(
    voidFractionModel,
    CGVoidFraction,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
CGVoidFraction::CGVoidFraction
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    voidFractionModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
	mesh_(particleCloud_.mesh()),
    verbose_(false),
    procBoundaryCorrection_(propsDict_.lookupOrDefault<Switch>("procBoundaryCorrection", false)),
    alphaMin_(readScalar(propsDict_.lookup("alphaMin"))),
	meshx_(readScalar(propsDict_.lookup("meshx"))),
	meshy_(readScalar(propsDict_.lookup("meshy"))),
	meshz_(readScalar(propsDict_.lookup("meshz"))),
	decx_(readScalar(propsDict_.lookup("decx"))),
	decy_(readScalar(propsDict_.lookup("decy"))),
	decz_(readScalar(propsDict_.lookup("decz"))),
	periodicx_(propsDict_.lookupOrDefault<Switch>("periodicx", false)),
    periodicy_(propsDict_.lookupOrDefault<Switch>("periodicy", false)),
    periodicz_(propsDict_.lookupOrDefault<Switch>("periodicz", false)),
	TopologyAccel_(propsDict_.lookupOrDefault<Switch>("TopologyAccel", false)),
    alphaLimited_(0),
    tooMuch_(0.0),
    interpolation_(false),
    cfdemUseOnly_(false)
{
    maxCellsPerParticle_ = 313;
    //particleCloud_.setMaxCellsPerParticle(29);

    if(alphaMin_ > 1 || alphaMin_ < 0.01){ FatalError<< "alphaMin should be < 1 and > 0.01 !!!" << abort(FatalError); }
    if (propsDict_.found("interpolation")){
        interpolation_=true;
        Warning << "interpolation for CGVoidFraction does not yet work correctly!" << endl;
        //Info << "Using interpolated voidfraction field - do not use this in combination with interpolation in drag model!"<< endl;
    }

    checkWeightNporosity(propsDict_);

    if (propsDict_.found("verbose")) verbose_=true;

    if (propsDict_.found("cfdemUseOnly"))
    {
        cfdemUseOnly_ = readBool(propsDict_.lookup("cfdemUseOnly"));
    }

    // check if settings are consistent with locate model selected
    if (procBoundaryCorrection_)
    {
        if(!(particleCloud_.locateM().type()=="engineIB"))
        {
            FatalError << typeName << ": You are requesting procBoundaryCorrection, this requires the use of engineIB!\n"
                       << abort(FatalError);
        }
    } else {
        if(particleCloud_.locateM().type()=="engineIB")
        {
            FatalError << typeName << ": You are using engineIB, this requires using procBoundaryCorrection=true!\n"
                       << abort(FatalError);
            //Warning << "You are trying to use engineIB, this requires using procBoundaryCorrection=true\n"
            //        << "  procBoundaryCorrection will be used!\n" << endl;
            //procBoundaryCorrection_ = true;
        }//
    }
	
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

CGVoidFraction::~CGVoidFraction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label CGVoidFraction::calculateLocalGridID(
    const vector& subPosition,
    const boundBox& globalBb,
    int meshx, int meshy, int meshz,
    scalar decx_, scalar decy_, scalar decz_) const
{
    // 网格尺寸计算
    scalar dx = (globalBb.max()[0] - globalBb.min()[0]) / meshx; // x 方向网格尺寸
    scalar dy = (globalBb.max()[1] - globalBb.min()[1]) / meshy; // y 方向网格尺寸
    scalar dz = (globalBb.max()[2] - globalBb.min()[2]) / meshz; // z 方向网格尺寸

    // 转换 decx_, decy_, decz_ 为整数
    int decx = static_cast<int>(decx_);
    int decy = static_cast<int>(decy_);
    int decz = static_cast<int>(decz_);

    // 当前处理器编号及在 x, y, z 三个方向的索引
    int procID = Pstream::myProcNo();
    int procX = procID % decx;
    int procY = (procID / decx) % decy;
    int procZ = procID / (decx * decy);

    // 局部处理器负责的网格数量
    int ProMeshNumX = meshx / decx;
    int ProMeshNumY = meshy / decy;
    int ProMeshNumZ = meshz / decz;

    // 局部处理器范围
    scalar localProcessorXMin = globalBb.min()[0] + procX * ProMeshNumX * dx;
    scalar localProcessorXMax = localProcessorXMin + ProMeshNumX * dx;

    scalar localProcessorYMin = globalBb.min()[1] + procY * ProMeshNumY * dy;
    scalar localProcessorYMax = localProcessorYMin + ProMeshNumY * dy;

    scalar localProcessorZMin = globalBb.min()[2] + procZ * ProMeshNumZ * dz;
    scalar localProcessorZMax = localProcessorZMin + ProMeshNumZ * dz;

    // 检查点是否在当前处理器范围内
    if (subPosition[0] < localProcessorXMin || subPosition[0] >= localProcessorXMax ||
        subPosition[1] < localProcessorYMin || subPosition[1] >= localProcessorYMax ||
        subPosition[2] < localProcessorZMin || subPosition[2] >= localProcessorZMax)
    {
        return -1; // 点不在当前处理器范围
    }

    // 计算局部网格索引
    int I = static_cast<int>((subPosition[0] - localProcessorXMin) / dx);
    int J = static_cast<int>((subPosition[1] - localProcessorYMin) / dy);
    int K = static_cast<int>((subPosition[2] - localProcessorZMin) / dz);

    // 根据索引公式计算局部网格 ID
    //return K * (ProMeshNumX * ProMeshNumY) + (ProMeshNumX - 1 - I) * ProMeshNumY + J;
	return K * (ProMeshNumX * ProMeshNumY) + J * ProMeshNumX + I;
}



void CGVoidFraction::setvoidFraction(double** const& mask,double**& voidfractions,double**& particleWeights,double**& particleVolumes, double**& particleV) const
{
    if(cfdemUseOnly_)
        reAllocArrays(particleCloud_.numberOfParticles());
    else
        reAllocArrays();
    voidfractionNext_ == dimensionedScalar("one", voidfractionNext_.dimensions(), 1.);//0826
	//for(int aa=0; aa< 1000; aa++)
	//{Info << " aa=" << aa << endl; Info << " voidfractionNext_=" << voidfractionNext_[aa] << endl; }
	
    vector position(0.,0.,0.);
    label cellID = -1;
    scalar radius(-1.);
    scalar volume(0.);
    scalar cellVol(0.);
    scalar scaleVol= weight();
    scalar scaleRadius = cbrt(porosity());
	scalar ds(0);
    const boundBox& globalBb = particleCloud_.mesh().bounds();
    label partCellId = -1 ;	
////////////////////////////////////布置卫星点系数//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	scalar r[] = { 0.5,1,1.5,2,2.5,3,3.5,4 };//定义卫星点层数 
	vector offsets[numberOfMarkerPoints];// add
    int m = 0;
    offsets[m][0] = offsets[m][1] = offsets[m][2] = 0.0;//初始化，并给圆心赋值
    m ++;//m=卫星点	
	double phipi=(sqrt(5.)-1)/2;
	    //第0.5层6
        for (int j = 1; j <= 6; j++)
	    {
			offsets[m][2] = r[0]*((2*static_cast<double>(j)-1)/6-1);
			offsets[m][0] = r[0]*sqrt(1-((2*static_cast<double>(j)-1)/6-1)*((2*static_cast<double>(j)-1)/6-1))*Foam::cos(2*M_PI*phipi*static_cast<double>(j));
            offsets[m][1] = r[0]*sqrt(1-((2*static_cast<double>(j)-1)/6-1)*((2*static_cast<double>(j)-1)/6-1))*Foam::sin(2*M_PI*phipi*static_cast<double>(j));
			m++;
		}
        //第一层6
        for (int j = 1; j <= 6; j++)
	    {
			offsets[m][2] = r[1]*((2*static_cast<double>(j)-1)/6-1);
			offsets[m][0] = r[1]*sqrt(1-((2*static_cast<double>(j)-1)/6-1)*((2*static_cast<double>(j)-1)/6-1))*Foam::cos(2*M_PI*phipi*static_cast<double>(j));
            offsets[m][1] = r[1]*sqrt(1-((2*static_cast<double>(j)-1)/6-1)*((2*static_cast<double>(j)-1)/6-1))*Foam::sin(2*M_PI*phipi*static_cast<double>(j));
			m++;
		}
		//第1.5层14
        for (int j = 1; j <= 14; j++)
	    {
			offsets[m][2] = r[2]*((2*static_cast<double>(j)-1)/14-1);
			offsets[m][0] = r[2]*sqrt(1-((2*static_cast<double>(j)-1)/14-1)*((2*static_cast<double>(j)-1)/14-1))*Foam::cos(2*M_PI*phipi*static_cast<double>(j));
            offsets[m][1] = r[2]*sqrt(1-((2*static_cast<double>(j)-1)/14-1)*((2*static_cast<double>(j)-1)/14-1))*Foam::sin(2*M_PI*phipi*static_cast<double>(j));
			m++;
		}
		//第二层24
        for (int j = 1; j <= 24; j++) 
        {
			offsets[m][2] = r[3]*((2*static_cast<double>(j)-1)/24-1);
			offsets[m][0] = r[3]*sqrt(1-((2*static_cast<double>(j)-1)/24-1)*((2*static_cast<double>(j)-1)/24-1))*Foam::cos(2*M_PI*phipi*static_cast<double>(j));
            offsets[m][1] = r[3]*sqrt(1-((2*static_cast<double>(j)-1)/24-1)*((2*static_cast<double>(j)-1)/24-1))*Foam::sin(2*M_PI*phipi*static_cast<double>(j));
            m ++;//Pout << " m2=" << m << endl;
        }
		//第2.5层38
        for (int j = 1; j <= 38; j++) 
        {
			offsets[m][2] = r[4]*((2*static_cast<double>(j)-1)/38-1);
			offsets[m][0] = r[4]*sqrt(1-((2*static_cast<double>(j)-1)/38-1)*((2*static_cast<double>(j)-1)/38-1))*Foam::cos(2*M_PI*phipi*static_cast<double>(j));
            offsets[m][1] = r[4]*sqrt(1-((2*static_cast<double>(j)-1)/38-1)*((2*static_cast<double>(j)-1)/38-1))*Foam::sin(2*M_PI*phipi*static_cast<double>(j));
            m ++;//Pout << " m2=" << m << endl;
        }
		//第三层54
		for (int j = 1; j <= 54; j++) 
        {
			offsets[m][2] = r[5]*((2*static_cast<double>(j)-1)/54-1);
			offsets[m][0] = r[5]*sqrt(1-((2*static_cast<double>(j)-1)/54-1)*((2*static_cast<double>(j)-1)/54-1))*Foam::cos(2*M_PI*phipi*static_cast<double>(j));
            offsets[m][1] = r[5]*sqrt(1-((2*static_cast<double>(j)-1)/54-1)*((2*static_cast<double>(j)-1)/54-1))*Foam::sin(2*M_PI*phipi*static_cast<double>(j));
            m ++;//Pout << " m3=" << m << endl;
        }
		//第三层74
		for (int j = 1; j <= 74; j++) 
        {
			offsets[m][2] = r[6]*((2*static_cast<double>(j)-1)/74-1);
			offsets[m][0] = r[6]*sqrt(1-((2*static_cast<double>(j)-1)/74-1)*((2*static_cast<double>(j)-1)/74-1))*Foam::cos(2*M_PI*phipi*static_cast<double>(j));
            offsets[m][1] = r[6]*sqrt(1-((2*static_cast<double>(j)-1)/74-1)*((2*static_cast<double>(j)-1)/74-1))*Foam::sin(2*M_PI*phipi*static_cast<double>(j));
            m ++;//Pout << " m3=" << m << endl;
        }
		//第四层96
		for (int j = 1; j <= 96; j++) 
        {
			offsets[m][2] = r[7]*((2*static_cast<double>(j)-1)/96-1);
			offsets[m][0] = r[7]*sqrt(1-((2*static_cast<double>(j)-1)/96-1)*((2*static_cast<double>(j)-1)/96-1))*Foam::cos(2*M_PI*phipi*static_cast<double>(j));
            offsets[m][1] = r[7]*sqrt(1-((2*static_cast<double>(j)-1)/96-1)*((2*static_cast<double>(j)-1)/96-1))*Foam::sin(2*M_PI*phipi*static_cast<double>(j));
            m ++;//Pout << " m4=" << m << endl;
        }
//#pragma omp parallel for
if (TopologyAccel_)
{
	for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {

        if(!checkParticleType(index)) continue; //skip this particle if not correct type

        //if(mask[index][0])
        //{
            // reset

            for(int subcell=0;subcell<cellsPerParticle_[index][0];subcell++)
            {
                particleWeights[index][subcell] = 0.;
                particleVolumes[index][subcell] = 0.;
            }
            particleV[index][0] = 0.;

            cellsPerParticle_[index][0] = 1.;//cellsperparticle是颗粒占据（重叠）的网格数
            position = particleCloud_.position(index);
            cellID = particleCloud_.cellIDs()[index][0];
			radius = particleRadius(index);//particleCloud_.radius(index);
            volume = Vp(index,radius,scaleVol);
            radius *= scaleRadius;
            cellVol = 0.;
            ds = 2*particleCloud_.radius(index);
            //--variables for sub-search
            int nPoints = numberOfMarkerPoints; //numberOfMarkerPoints 头文件定义
            int nNotFound=0,nUnEqual=0,nTotal=0;
            vector offset(0.,0.,0.);
            int cellsSet = 0;
            label cellWithCenter(-1);
			//在OpenFOAM中，label 是一种特定的整数类型，用于标识单元、面、边、点等网格元素的标签
            //label 类型在OpenFOAM中的定义通常是一个32位整数（即int类型），但这可能取决于特定版本或编译选项
			//label 类型的变量是OpenFOAM中用于标识网格元素的一种特殊整数类型，具有上述的特点和用途
			if(procBoundaryCorrection_)
            {
                // switch off cellIDs for force calc if steming from parallel search success
                cellWithCenter = particleCloud_.locateM().findSingleCell(position,cellID);
                particleCloud_.cellIDs()[index][0] = cellWithCenter;
            }
	labelList transffi(313,0);
	int count=0;
	int cc=0;
/*///////////////局部-全局网格ID/////////////////////////////////////////////////////
#include "globalMeshData.H"
#include "globalIndex.H"
const polyMesh& mesh = particleCloud_.mesh();
Foam::globalMeshData globalData(mesh); // 创建 globalMeshData 对象
const Foam::globalIndex& globalCellIndex = globalData.globalBoundaryCellNumbering(); // 获取全局网格编号工具
label globalID = globalCellIndex.toGlobal(cellID);
// 输出结果
Pout << "Processor " << Foam::Pstream::myProcNo()
     << ", Local Cell ID: " << cellID
     << ", Global Cell ID: " << globalID << Foam::endl;
	 	
#include "Pstream.H"
#include "globalMeshData.H"
#include "globalIndex.H"
#include "HashSet.H"

const polyMesh& mesh = particleCloud_.mesh();
const Foam::globalMeshData globalData(mesh);

// 获取全局编号工具
const Foam::globalIndex& globalCellIndex = globalData.globalBoundaryCellNumbering();

// 获取偏移量
const Foam::labelList& ooffsets = const_cast<Foam::globalIndex&>(globalCellIndex).offsets();

Foam::Pout << "Processor " << Foam::Pstream::myProcNo()
           << ", Offsets: " << ooffsets << Foam::endl;


/*for (Foam::label localID = 0; localID < mesh.nCells(); ++localID)
{
    Foam::label globalID = globalCellIndex.toGlobal(localID);

    // 判断当前处理器是否是主处理器
    if (globalCellIndex.whichProcID(globalID) == Foam::Pstream::myProcNo())
    {
        Foam::Pout << "Processor " << Foam::Pstream::myProcNo()
                   << ", Local Cell ID: " << localID
                   << ", Global Cell ID: " << globalID << " (Master)" << Foam::endl;
    }
}*/


/*const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();
forAll(boundaryMesh, patchi)
{
    const polyPatch& patch = boundaryMesh[patchi];

    if (isA<processorPolyPatch>(patch))
    {
        const processorPolyPatch& procPatch = refCast<const processorPolyPatch>(patch);
        const labelList& faceCells = procPatch.faceCells();

        // 遍历共享单元
        forAll(faceCells, i)
        {
            Foam::Pout << "Shared Cell ID: " << faceCells[i] << Foam::endl;
        }
    }
}*/

//////////////////////////////////////////////////////////////////
 if (cellID >= 0)  
 {
    cellVol = particleCloud_.mesh().V()[cellID];
	scalar lostweight=0;
  for(int i = 0; i < numberOfMarkerPoints; i++) //遍历所有卫星点，计算所有卫星点的偏移量
  {

//////////////////////////////////////////CGprocess//////////////////////////////////////////////////////////////////						                 
    // locate subPoint
	                vector offset(0.,0.,0.);
					offset = particleRadius(index)*offsets[i];
					scalar distance = mag (offset);
					vector subPosition = particleCloud_.position(index) + offset;
					cellVol = particleCloud_.mesh().V()[cellID];
					scalar wl = 0;
					scalar layerR=0;

		// 对 subPosition 的每个方向进行周期性修正，iDir: 0->x, 1->y, 2->z
		for (int iDir = 0; iDir < 3; iDir++)
		{
 		   // 判断当前方向是否启用周期性边界
 		   bool isPeriodic = false;
 		   if (iDir == 0 && periodicx_)
 		   {
 		       isPeriodic = true;
 		   }
 		   else if (iDir == 1 && periodicy_)
 		   {
 		       isPeriodic = true;
 		   }
 		   else if (iDir == 2 && periodicz_)
 		   {
  		      isPeriodic = true;
 		   }
		
  		    if (isPeriodic)
  		   {
  		      if (subPosition[iDir] > globalBb.max()[iDir])
   		     {
           		 subPosition[iDir] -= (globalBb.max()[iDir] - globalBb.min()[iDir]);
        		}
 		       else if (subPosition[iDir] < globalBb.min()[iDir])
   		     {
   		         subPosition[iDir] += (globalBb.max()[iDir] - globalBb.min()[iDir]);
  		      }
   		    }
		}

	  
					
                        if (i == 0)
						{ wl=0.130856838989584*exp(-distance*distance/(8*ds*ds)); }
					
					else if (i >= 1 && i <= 6)
					{ wl=0.130856838989584/6*exp(-distance*distance/(8*ds*ds)); layerR = 0.5 ;}
				
					    else if (i >= 7 && i <= 12)
						{ wl=0.130856838989584/6*exp(-distance*distance/(8*ds*ds)); layerR = 1 ;}
					
					else if (i >= 13 && i <= 26)
					{ wl=0.130856838989584/14*exp(-distance*distance/(8*ds*ds)); layerR = 1.5 ;}
				
					    else if (i >= 27 && i <= 50)
						{ wl=0.130856838989584/24*exp(-distance*distance/(8*ds*ds)); layerR = 2 ;}
					
					else if (i >= 51 && i <= 88)
					{ wl=0.130856838989584/38*exp(-distance*distance/(8*ds*ds)); layerR = 2.5 ;}
				
					    else if (i >= 89 && i <= 142)
						{ wl=0.130856838989584/54*exp(-distance*distance/(8*ds*ds)); layerR = 3 ;}
					
					else if (i >= 143 && i <= 216)
					{ wl=0.130856838989584/74*exp(-distance*distance/(8*ds*ds)); layerR = 3.5 ;}
				
					    else if (i >= 217 && i <= 312)
						{ wl=0.130856838989584/96*exp(-distance*distance/(8*ds*ds)); layerR = 4 ;}
	//if (TopologyAccel_)			
    label partCellId = calculateLocalGridID(subPosition, globalBb, meshx_, meshy_, meshz_, decx_, decy_, decz_); 
    //else
	//{label partCellId = particleCloud_.locateM().findSingleCell(subPosition,cellID); }	
	//label partCellId = particleCloud_.locateM().findSingleCell(subPosition,cellID);
	//label partCellIdtrue = particleCloud_.locateM().findSingleCell(subPosition,cellID);
    //Pout << "partCellId true =" << partCellIdtrue << " ----partCellId TST = " << partCellId <<  endl;
     if (verbose_ && index==0)  	 
		 {
				label partCellIdtrue = particleCloud_.locateM().findSingleCell(subPosition,cellID);
                Pout << "partCellId true =" << partCellIdtrue << "--V-S--partCellId TST = " << partCellId <<  endl;
		 }	
////////////////////////////virtual subparticle method liuyuxiang2024.9.16	
    if (partCellId >= 0)  // subPoint is in domain 检查卫星点是不是在计算域内
    {
        // update voidfraction for each particle read
        scalar partCellVol = particleCloud_.mesh().V()[partCellId];//提取这个包含卫星点的单元的体积
        scalar particleVolume = volume*wl;//卫星点加权体积
      
        scalar newAlpha = voidfractionNext_[partCellId]- particleVolume / partCellVol;
		
        if(newAlpha > alphaMin_) voidfractionNext_[partCellId] = newAlpha;
        else
        {
            voidfractionNext_[partCellId] = alphaMin_;
            tooMuch_ += (alphaMin_-newAlpha) * partCellVol;
        }
        cellsSet++;   
//--------set sub weight--------////--------set sub weight--------////--------set sub weight--------////--------set sub weight--------//
        bool createNew = true;
        label storeInIndex=0;
        for(int i=0; i < cellsPerParticle_[index][0] ; i++)
        {
            if(partCellId == particleCloud_.cellIDs()[index][i]) 
            {
                storeInIndex = i;
                createNew = false;
                break;
            }
        }
        if(createNew)
        {
            cellsPerParticle_[index][0] ++;
            storeInIndex = cellsPerParticle_[index][0]-1;
            particleCloud_.cellIDs()[index][storeInIndex] = partCellId;
        }
        particleWeights[index][storeInIndex] += wl; 
        particleVolumes[index][storeInIndex] += particleVolume;  
        particleV[index][0] += particleVolume;	
    }
	if (partCellId < 0) 
	{    
        ////////////////////////////////////////first physical boundary
		if (subPosition[0] > globalBb.max()[0] || subPosition[1] > globalBb.max()[1] || subPosition[2] > globalBb.max()[2] || subPosition[0] < globalBb.min()[0] || subPosition[1] < globalBb.min()[1] || subPosition[2] < globalBb.min()[2])
	    {    
            while (partCellId < 0)//&& iterations < iterationslimit)	//virtual sub-point
            { 
				if (layerR <= 0) //treatment for point inside processor boundary!	(very little)
				{
					break;
				}		
	            for (int ivsp = 0 ; ivsp < 3 ; ivsp++)
	            {
                    offset[ivsp] = offset[ivsp]*(layerR-0.5)/(layerR);// layerR must not be 0
	            }	
	            layerR=layerR-0.5;
	            subPosition = particleCloud_.position(index) + offset;	
	            //partCellId = particleCloud_.locateM().findSingleCell(subPosition,cellID);
				//partCellId = calculateLocalGridID(subPosition, globalBb, meshx_, meshy_, meshz_, decx_, decy_, decz_);
    			//if (TopologyAccel_)			
   			    partCellId = calculateLocalGridID(subPosition, globalBb, meshx_, meshy_, meshz_, decx_, decy_, decz_); 
   			    //else
				//{partCellId = particleCloud_.locateM().findSingleCell(subPosition,cellID); }	
            }
			
			if (partCellId >= 0) //////// move the virtual point into computational domain
			{   
			    scalar partCellVol = particleCloud_.mesh().V()[partCellId];//提取这个包含卫星点的单元的体积
                scalar particleVolume = volume*wl;//卫星点加权体积
                scalar newAlpha = voidfractionNext_[partCellId]- particleVolume / partCellVol;
		        
        	    if(newAlpha > alphaMin_) voidfractionNext_[partCellId] = newAlpha;
        	    else
        	    {
        	        voidfractionNext_[partCellId] = alphaMin_;
        	        tooMuch_ += (alphaMin_-newAlpha) * partCellVol;
        	    }
        	    cellsSet++; 
//--------set sub weight--------////--------set sub weight--------////--------set sub weight--------////--------set sub weight--------//
				bool createNew = true;
 				label storeInIndex=0;
  				for(int i=0; i < cellsPerParticle_[index][0] ; i++)
  				{
    				if(partCellId == particleCloud_.cellIDs()[index][i]) 
    				{
      				    storeInIndex = i;
     				    createNew = false;
      				    break;
     				}
     			}
     			if(createNew)
     			{
      				cellsPerParticle_[index][0] ++;
      				storeInIndex = cellsPerParticle_[index][0]-1;
      				particleCloud_.cellIDs()[index][storeInIndex] = partCellId;
      			}
      			particleWeights[index][storeInIndex] += wl; 
      			particleVolumes[index][storeInIndex] += particleVolume;  
       			particleV[index][0] += particleVolume;
			}
	    }//end physical boundary revision (out of physcial boundary!)
		else // inside the physcial boundary but not on this processor!
		{ 
		    //权重修正，体积分数用reduce准确传递，权重只在该处理器修正
			lostweight= lostweight+wl; //不修正体积分数只修正权重，因为体积分数在后面的reduce处理！！！！！！
		    transffi[cc]=i;
		    cc=cc+1;
		    count=count+1;
		}
    }//end for not main processor (partcellid not found!)
 }// end subpoint loop
    for(int subcellss=0;subcellss<cellsPerParticle_[index][0];subcellss++)
    {
	  particleWeights[index][subcellss]=particleWeights[index][subcellss]/(1-lostweight);
	  particleVolumes[index][subcellss]=particleVolumes[index][subcellss]/(1-lostweight);
    }
	particleV[index][0]=particleV[index][0]/(1-lostweight);
}//end if in cell main processor

	 reduce(transffi, maxOp<List<label>>());// test OK transffi=没找到的partID对应的i（点云次序）
	 reduce(count, maxOp<int>());// test OK transffi=没找到的partID对应的i（点云次序）
if (cellID < 0)
{
	for (cc=0;cc<count;cc++) //在非主处理器上遍历不在主处理器上的卫星点
	{
		int i=transffi[cc]; // transffi存的是i（从小到大），点云次序==等价于主处理器上的i 只是这里的i是没被搜索到的
		vector offset(0.,0.,0.);
		offset = particleRadius(index)*offsets[i];
		scalar distance = mag (offset);
		vector subPosition = particleCloud_.position(index) + offset;
		cellVol = particleCloud_.mesh().V()[cellID];
		scalar wl = 0;
		scalar layerR=0;
// 对 subPosition 的每个方向进行周期性修正，iDir: 0->x, 1->y, 2->z
		for (int iDir = 0; iDir < 3; iDir++)
		{
 		   // 判断当前方向是否启用周期性边界
 		   bool isPeriodic = false;
 		   if (iDir == 0 && periodicx_)
 		   {
 		       isPeriodic = true;
 		   }
 		   else if (iDir == 1 && periodicy_)
 		   {
 		       isPeriodic = true;
 		   }
 		   else if (iDir == 2 && periodicz_)
 		   {
  		      isPeriodic = true;
 		   }
		
 		   // 如果该方向为周期性，则进行坐标“包裹”处理
  		    if (isPeriodic)
  		   {
  		      if (subPosition[iDir] > globalBb.max()[iDir])
   		     {
           		 subPosition[iDir] -= (globalBb.max()[iDir] - globalBb.min()[iDir]);
        		}
 		       else if (subPosition[iDir] < globalBb.min()[iDir])
   		     {
   		         subPosition[iDir] += (globalBb.max()[iDir] - globalBb.min()[iDir]);
  		      }
   		    }
		}
                        if (i == 0)
						{ wl=0.130856838989584*exp(-distance*distance/(8*ds*ds)); }
					
					else if (i >= 1 && i <= 6)
					{ wl=0.130856838989584/6*exp(-distance*distance/(8*ds*ds)); layerR = 0.5 ;}
				
					    else if (i >= 7 && i <= 12)
						{ wl=0.130856838989584/6*exp(-distance*distance/(8*ds*ds)); layerR = 1 ;}
					
					else if (i >= 13 && i <= 26)
					{ wl=0.130856838989584/14*exp(-distance*distance/(8*ds*ds)); layerR = 1.5 ;}
				
					    else if (i >= 27 && i <= 50)
						{ wl=0.130856838989584/24*exp(-distance*distance/(8*ds*ds)); layerR = 2 ;}
					
					else if (i >= 51 && i <= 88)
					{ wl=0.130856838989584/38*exp(-distance*distance/(8*ds*ds)); layerR = 2.5 ;}
				
					    else if (i >= 89 && i <= 142)
						{ wl=0.130856838989584/54*exp(-distance*distance/(8*ds*ds)); layerR = 3 ;}
					
					else if (i >= 143 && i <= 216)
					{ wl=0.130856838989584/74*exp(-distance*distance/(8*ds*ds)); layerR = 3.5 ;}
				
					    else if (i >= 217 && i <= 312)
						{ wl=0.130856838989584/96*exp(-distance*distance/(8*ds*ds)); layerR = 4 ;}
                    //label partCellId = particleCloud_.locateM().findSingleCell(subPosition,cellID);
                    //label partCellId = calculateLocalGridID(subPosition, globalBb, meshx_, meshy_, meshz_, decx_, decy_, decz_);
				//if (TopologyAccel_)			
   			    label partCellId = calculateLocalGridID(subPosition, globalBb, meshx_, meshy_, meshz_, decx_, decy_, decz_); 
   			    //else
				//{label partCellId = particleCloud_.locateM().findSingleCell(subPosition,cellID); }
	  if (partCellId >= 0)///////////////非主处理器上的卫星点
	  {
        scalar partCellVol = particleCloud_.mesh().V()[partCellId];
        scalar particleVolume = volume*wl;
        scalar newAlpha = voidfractionNext_[partCellId]- particleVolume / partCellVol;
        if(newAlpha > alphaMin_) voidfractionNext_[partCellId] = newAlpha;
        else
        {
            voidfractionNext_[partCellId] = alphaMin_;
            tooMuch_ += (alphaMin_-newAlpha) * partCellVol;
        }
        cellsSet++; //cancel weighting! because the error of array only revise the voidfraction!
	  } 
	}// end subpoint非主处理器卫星点遍历		 
}// end 非主处理器	

	            if (verbose_ && index == 0)
				{
				   scalar sum=0; 
				   for (int lyx = 0; lyx < particleCloud_.cellsPerParticle()[index][0]; lyx++) //kernalcell是颗粒CG域占据的网格次序，storesubcellkernal才是实际的网格id
				   {
				       sum=sum+particleCloud_.particleWeights()[index][lyx]; //debug for mass conservation
                   }
				   Pout << "volume=" << volume << endl; 
				   Pout << "particleWeights[index][0]=" << particleWeights[index][0] << endl;
				   Pout << "particleVolumes[index][0]=" << particleVolumes[index][0] << endl;
				   Pout << "particleV[index][0]=" << particleV[index][0] << endl;
				   Pout << "sum=" << sum << endl;
				}

    }//// end loop all particles
}//// end TopologyAccel_
else
{		
    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {

        if(!checkParticleType(index)) continue; //skip this particle if not correct type

        //if(mask[index][0])
        //{
            // reset

            for(int subcell=0;subcell<cellsPerParticle_[index][0];subcell++)
            {
                particleWeights[index][subcell] = 0.;
                particleVolumes[index][subcell] = 0.;
            }
            particleV[index][0] = 0.;

            cellsPerParticle_[index][0] = 1.;//cellsperparticle是颗粒占据（重叠）的网格数
            position = particleCloud_.position(index);
            cellID = particleCloud_.cellIDs()[index][0];
			radius = particleRadius(index);//particleCloud_.radius(index);
            volume = Vp(index,radius,scaleVol);
            radius *= scaleRadius;
            cellVol = 0.;
            ds = 2*particleCloud_.radius(index);
            //--variables for sub-search
            int nPoints = numberOfMarkerPoints; //numberOfMarkerPoints 头文件定义
            int nNotFound=0,nUnEqual=0,nTotal=0;
            vector offset(0.,0.,0.);
            int cellsSet = 0;
            label cellWithCenter(-1);
			//在OpenFOAM中，label 是一种特定的整数类型，用于标识单元、面、边、点等网格元素的标签
            //label 类型在OpenFOAM中的定义通常是一个32位整数（即int类型），但这可能取决于特定版本或编译选项
			//label 类型的变量是OpenFOAM中用于标识网格元素的一种特殊整数类型，具有上述的特点和用途
			if(procBoundaryCorrection_)
            {
                // switch off cellIDs for force calc if steming from parallel search success
                cellWithCenter = particleCloud_.locateM().findSingleCell(position,cellID);
                particleCloud_.cellIDs()[index][0] = cellWithCenter;
            }
	labelList transffi(313,0);
	int count=0;
	int cc=0;
/*///////////////局部-全局网格ID/////////////////////////////////////////////////////
#include "globalMeshData.H"
#include "globalIndex.H"
const polyMesh& mesh = particleCloud_.mesh();
Foam::globalMeshData globalData(mesh); // 创建 globalMeshData 对象
const Foam::globalIndex& globalCellIndex = globalData.globalBoundaryCellNumbering(); // 获取全局网格编号工具
label globalID = globalCellIndex.toGlobal(cellID);
// 输出结果
Pout << "Processor " << Foam::Pstream::myProcNo()
     << ", Local Cell ID: " << cellID
     << ", Global Cell ID: " << globalID << Foam::endl;
	 	
#include "Pstream.H"
#include "globalMeshData.H"
#include "globalIndex.H"
#include "HashSet.H"

const polyMesh& mesh = particleCloud_.mesh();
const Foam::globalMeshData globalData(mesh);

// 获取全局编号工具
const Foam::globalIndex& globalCellIndex = globalData.globalBoundaryCellNumbering();

// 获取偏移量
const Foam::labelList& ooffsets = const_cast<Foam::globalIndex&>(globalCellIndex).offsets();

Foam::Pout << "Processor " << Foam::Pstream::myProcNo()
           << ", Offsets: " << ooffsets << Foam::endl;


/*for (Foam::label localID = 0; localID < mesh.nCells(); ++localID)
{
    Foam::label globalID = globalCellIndex.toGlobal(localID);

    // 判断当前处理器是否是主处理器
    if (globalCellIndex.whichProcID(globalID) == Foam::Pstream::myProcNo())
    {
        Foam::Pout << "Processor " << Foam::Pstream::myProcNo()
                   << ", Local Cell ID: " << localID
                   << ", Global Cell ID: " << globalID << " (Master)" << Foam::endl;
    }
}*/


/*const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();
forAll(boundaryMesh, patchi)
{
    const polyPatch& patch = boundaryMesh[patchi];

    if (isA<processorPolyPatch>(patch))
    {
        const processorPolyPatch& procPatch = refCast<const processorPolyPatch>(patch);
        const labelList& faceCells = procPatch.faceCells();

        // 遍历共享单元
        forAll(faceCells, i)
        {
            Foam::Pout << "Shared Cell ID: " << faceCells[i] << Foam::endl;
        }
    }
}*/

//////////////////////////////////////////////////////////////////
 if (cellID >= 0)  
 {
    cellVol = particleCloud_.mesh().V()[cellID];
	scalar lostweight=0;
  for(int i = 0; i < numberOfMarkerPoints; i++) //遍历所有卫星点，计算所有卫星点的偏移量
  {

//////////////////////////////////////////CGprocess//////////////////////////////////////////////////////////////////						                 
    // locate subPoint
	                vector offset(0.,0.,0.);
					offset = particleRadius(index)*offsets[i];
					scalar distance = mag (offset);
					vector subPosition = particleCloud_.position(index) + offset;
					cellVol = particleCloud_.mesh().V()[cellID];
					scalar wl = 0;
					scalar layerR=0;

		// 对 subPosition 的每个方向进行周期性修正，iDir: 0->x, 1->y, 2->z
		for (int iDir = 0; iDir < 3; iDir++)
		{
 		   // 判断当前方向是否启用周期性边界
 		   bool isPeriodic = false;
 		   if (iDir == 0 && periodicx_)
 		   {
 		       isPeriodic = true;
 		   }
 		   else if (iDir == 1 && periodicy_)
 		   {
 		       isPeriodic = true;
 		   }
 		   else if (iDir == 2 && periodicz_)
 		   {
  		      isPeriodic = true;
 		   }
		
  		    if (isPeriodic)
  		   {
  		      if (subPosition[iDir] > globalBb.max()[iDir])
   		     {
           		 subPosition[iDir] -= (globalBb.max()[iDir] - globalBb.min()[iDir]);
        		}
 		       else if (subPosition[iDir] < globalBb.min()[iDir])
   		     {
   		         subPosition[iDir] += (globalBb.max()[iDir] - globalBb.min()[iDir]);
  		      }
   		    }
		}

	  
					
                        if (i == 0)
						{ wl=0.130856838989584*exp(-distance*distance/(8*ds*ds)); }
					
					else if (i >= 1 && i <= 6)
					{ wl=0.130856838989584/6*exp(-distance*distance/(8*ds*ds)); layerR = 0.5 ;}
				
					    else if (i >= 7 && i <= 12)
						{ wl=0.130856838989584/6*exp(-distance*distance/(8*ds*ds)); layerR = 1 ;}
					
					else if (i >= 13 && i <= 26)
					{ wl=0.130856838989584/14*exp(-distance*distance/(8*ds*ds)); layerR = 1.5 ;}
				
					    else if (i >= 27 && i <= 50)
						{ wl=0.130856838989584/24*exp(-distance*distance/(8*ds*ds)); layerR = 2 ;}
					
					else if (i >= 51 && i <= 88)
					{ wl=0.130856838989584/38*exp(-distance*distance/(8*ds*ds)); layerR = 2.5 ;}
				
					    else if (i >= 89 && i <= 142)
						{ wl=0.130856838989584/54*exp(-distance*distance/(8*ds*ds)); layerR = 3 ;}
					
					else if (i >= 143 && i <= 216)
					{ wl=0.130856838989584/74*exp(-distance*distance/(8*ds*ds)); layerR = 3.5 ;}
				
					    else if (i >= 217 && i <= 312)
						{ wl=0.130856838989584/96*exp(-distance*distance/(8*ds*ds)); layerR = 4 ;}
	//if (TopologyAccel_)			
    //{label partCellId = calculateLocalGridID(subPosition, globalBb, meshx_, meshy_, meshz_, decx_, decy_, decz_); }
    //else
	label partCellId = particleCloud_.locateM().findSingleCell(subPosition,cellID); 	
	//label partCellId = particleCloud_.locateM().findSingleCell(subPosition,cellID);
	//label partCellIdtrue = particleCloud_.locateM().findSingleCell(subPosition,cellID);
    //Pout << "partCellId true =" << partCellIdtrue << " ----partCellId TST = " << partCellId <<  endl;
     if (verbose_ && index==0)  	 
		 {
				label partCellIdtrue = particleCloud_.locateM().findSingleCell(subPosition,cellID);
                Pout << "partCellId true =" << partCellIdtrue << "--V-S--partCellId TST = " << partCellId <<  endl;
		 }	
////////////////////////////virtual subparticle method liuyuxiang2024.9.16	
    if (partCellId >= 0)  // subPoint is in domain 检查卫星点是不是在计算域内
    {
        // update voidfraction for each particle read
        scalar partCellVol = particleCloud_.mesh().V()[partCellId];//提取这个包含卫星点的单元的体积
        scalar particleVolume = volume*wl;//卫星点加权体积
      
        scalar newAlpha = voidfractionNext_[partCellId]- particleVolume / partCellVol;
		
        if(newAlpha > alphaMin_) voidfractionNext_[partCellId] = newAlpha;
        else
        {
            voidfractionNext_[partCellId] = alphaMin_;
            tooMuch_ += (alphaMin_-newAlpha) * partCellVol;
        }
        cellsSet++;   
//--------set sub weight--------////--------set sub weight--------////--------set sub weight--------////--------set sub weight--------//
        bool createNew = true;
        label storeInIndex=0;
        for(int i=0; i < cellsPerParticle_[index][0] ; i++)
        {
            if(partCellId == particleCloud_.cellIDs()[index][i]) 
            {
                storeInIndex = i;
                createNew = false;
                break;
            }
        }
        if(createNew)
        {
            cellsPerParticle_[index][0] ++;
            storeInIndex = cellsPerParticle_[index][0]-1;
            particleCloud_.cellIDs()[index][storeInIndex] = partCellId;
        }
        particleWeights[index][storeInIndex] += wl; 
        particleVolumes[index][storeInIndex] += particleVolume;  
        particleV[index][0] += particleVolume;	
    }
	if (partCellId < 0) 
	{    
        ////////////////////////////////////////first physical boundary
		if (subPosition[0] > globalBb.max()[0] || subPosition[1] > globalBb.max()[1] || subPosition[2] > globalBb.max()[2] || subPosition[0] < globalBb.min()[0] || subPosition[1] < globalBb.min()[1] || subPosition[2] < globalBb.min()[2])
	    {    
            while (partCellId < 0)//&& iterations < iterationslimit)	//virtual sub-point
            { 
				if (layerR <= 0) //treatment for point inside processor boundary!	(very little)
				{
					break;
				}		
	            for (int ivsp = 0 ; ivsp < 3 ; ivsp++)
	            {
                    offset[ivsp] = offset[ivsp]*(layerR-0.5)/(layerR);// layerR must not be 0
	            }	
	            layerR=layerR-0.5;
	            subPosition = particleCloud_.position(index) + offset;	
	            //partCellId = particleCloud_.locateM().findSingleCell(subPosition,cellID);
				//partCellId = calculateLocalGridID(subPosition, globalBb, meshx_, meshy_, meshz_, decx_, decy_, decz_);
    			//if (TopologyAccel_)			
   			    //{partCellId = calculateLocalGridID(subPosition, globalBb, meshx_, meshy_, meshz_, decx_, decy_, decz_); }
   			    //else
				{partCellId = particleCloud_.locateM().findSingleCell(subPosition,cellID); }	
            }
			
			if (partCellId >= 0) //////// move the virtual point into computational domain
			{   
			    scalar partCellVol = particleCloud_.mesh().V()[partCellId];//提取这个包含卫星点的单元的体积
                scalar particleVolume = volume*wl;//卫星点加权体积
                scalar newAlpha = voidfractionNext_[partCellId]- particleVolume / partCellVol;
		        
        	    if(newAlpha > alphaMin_) voidfractionNext_[partCellId] = newAlpha;
        	    else
        	    {
        	        voidfractionNext_[partCellId] = alphaMin_;
        	        tooMuch_ += (alphaMin_-newAlpha) * partCellVol;
        	    }
        	    cellsSet++; 
//--------set sub weight--------////--------set sub weight--------////--------set sub weight--------////--------set sub weight--------//
				bool createNew = true;
 				label storeInIndex=0;
  				for(int i=0; i < cellsPerParticle_[index][0] ; i++)
  				{
    				if(partCellId == particleCloud_.cellIDs()[index][i]) 
    				{
      				    storeInIndex = i;
     				    createNew = false;
      				    break;
     				}
     			}
     			if(createNew)
     			{
      				cellsPerParticle_[index][0] ++;
      				storeInIndex = cellsPerParticle_[index][0]-1;
      				particleCloud_.cellIDs()[index][storeInIndex] = partCellId;
      			}
      			particleWeights[index][storeInIndex] += wl; 
      			particleVolumes[index][storeInIndex] += particleVolume;  
       			particleV[index][0] += particleVolume;
			}
	    }//end physical boundary revision (out of physcial boundary!)
		else // inside the physcial boundary but not on this processor!
		{ 
		    //权重修正，体积分数用reduce准确传递，权重只在该处理器修正
			lostweight= lostweight+wl; //不修正体积分数只修正权重，因为体积分数在后面的reduce处理！！！！！！
		    transffi[cc]=i;
		    cc=cc+1;
		    count=count+1;
		}
    }//end for not main processor (partcellid not found!)
 }// end subpoint loop
    for(int subcellss=0;subcellss<cellsPerParticle_[index][0];subcellss++)
    {
	  particleWeights[index][subcellss]=particleWeights[index][subcellss]/(1-lostweight);
	  particleVolumes[index][subcellss]=particleVolumes[index][subcellss]/(1-lostweight);
    }
	particleV[index][0]=particleV[index][0]/(1-lostweight);
}//end if in cell main processor

	 reduce(transffi, maxOp<List<label>>());// test OK transffi=没找到的partID对应的i（点云次序）
	 reduce(count, maxOp<int>());// test OK transffi=没找到的partID对应的i（点云次序）
if (cellID < 0)
{
	for (cc=0;cc<count;cc++) //在非主处理器上遍历不在主处理器上的卫星点
	{
		int i=transffi[cc]; // transffi存的是i（从小到大），点云次序==等价于主处理器上的i 只是这里的i是没被搜索到的
		vector offset(0.,0.,0.);
		offset = particleRadius(index)*offsets[i];
		scalar distance = mag (offset);
		vector subPosition = particleCloud_.position(index) + offset;
		cellVol = particleCloud_.mesh().V()[cellID];
		scalar wl = 0;
		scalar layerR=0;
// 对 subPosition 的每个方向进行周期性修正，iDir: 0->x, 1->y, 2->z
		for (int iDir = 0; iDir < 3; iDir++)
		{
 		   // 判断当前方向是否启用周期性边界
 		   bool isPeriodic = false;
 		   if (iDir == 0 && periodicx_)
 		   {
 		       isPeriodic = true;
 		   }
 		   else if (iDir == 1 && periodicy_)
 		   {
 		       isPeriodic = true;
 		   }
 		   else if (iDir == 2 && periodicz_)
 		   {
  		      isPeriodic = true;
 		   }
		
 		   // 如果该方向为周期性，则进行坐标“包裹”处理
  		    if (isPeriodic)
  		   {
  		      if (subPosition[iDir] > globalBb.max()[iDir])
   		     {
           		 subPosition[iDir] -= (globalBb.max()[iDir] - globalBb.min()[iDir]);
        		}
 		       else if (subPosition[iDir] < globalBb.min()[iDir])
   		     {
   		         subPosition[iDir] += (globalBb.max()[iDir] - globalBb.min()[iDir]);
  		      }
   		    }
		}
                        if (i == 0)
						{ wl=0.130856838989584*exp(-distance*distance/(8*ds*ds)); }
					
					else if (i >= 1 && i <= 6)
					{ wl=0.130856838989584/6*exp(-distance*distance/(8*ds*ds)); layerR = 0.5 ;}
				
					    else if (i >= 7 && i <= 12)
						{ wl=0.130856838989584/6*exp(-distance*distance/(8*ds*ds)); layerR = 1 ;}
					
					else if (i >= 13 && i <= 26)
					{ wl=0.130856838989584/14*exp(-distance*distance/(8*ds*ds)); layerR = 1.5 ;}
				
					    else if (i >= 27 && i <= 50)
						{ wl=0.130856838989584/24*exp(-distance*distance/(8*ds*ds)); layerR = 2 ;}
					
					else if (i >= 51 && i <= 88)
					{ wl=0.130856838989584/38*exp(-distance*distance/(8*ds*ds)); layerR = 2.5 ;}
				
					    else if (i >= 89 && i <= 142)
						{ wl=0.130856838989584/54*exp(-distance*distance/(8*ds*ds)); layerR = 3 ;}
					
					else if (i >= 143 && i <= 216)
					{ wl=0.130856838989584/74*exp(-distance*distance/(8*ds*ds)); layerR = 3.5 ;}
				
					    else if (i >= 217 && i <= 312)
						{ wl=0.130856838989584/96*exp(-distance*distance/(8*ds*ds)); layerR = 4 ;}
                    //label partCellId = particleCloud_.locateM().findSingleCell(subPosition,cellID);
                    //label partCellId = calculateLocalGridID(subPosition, globalBb, meshx_, meshy_, meshz_, decx_, decy_, decz_);
				//if (TopologyAccel_)			
   			    //{label partCellId = calculateLocalGridID(subPosition, globalBb, meshx_, meshy_, meshz_, decx_, decy_, decz_); }
   			    //else
				label partCellId = particleCloud_.locateM().findSingleCell(subPosition,cellID); 
	  if (partCellId >= 0)///////////////非主处理器上的卫星点
	  {
        scalar partCellVol = particleCloud_.mesh().V()[partCellId];
        scalar particleVolume = volume*wl;
        scalar newAlpha = voidfractionNext_[partCellId]- particleVolume / partCellVol;
        if(newAlpha > alphaMin_) voidfractionNext_[partCellId] = newAlpha;
        else
        {
            voidfractionNext_[partCellId] = alphaMin_;
            tooMuch_ += (alphaMin_-newAlpha) * partCellVol;
        }
        cellsSet++; //cancel weighting! because the error of array only revise the voidfraction!
	  } 
	}// end subpoint非主处理器卫星点遍历		 
}// end 非主处理器	

	            if (verbose_ && index == 0)
				{
				   scalar sum=0; 
				   for (int lyx = 0; lyx < particleCloud_.cellsPerParticle()[index][0]; lyx++) //kernalcell是颗粒CG域占据的网格次序，storesubcellkernal才是实际的网格id
				   {
				       sum=sum+particleCloud_.particleWeights()[index][lyx]; //debug for mass conservation
                   }
				   Pout << "volume=" << volume << endl; 
				   Pout << "particleWeights[index][0]=" << particleWeights[index][0] << endl;
				   Pout << "particleVolumes[index][0]=" << particleVolumes[index][0] << endl;
				   Pout << "particleV[index][0]=" << particleV[index][0] << endl;
				   Pout << "sum=" << sum << endl;
				}

    }//// end loop all particles
}





///////////////new loop test/////////////////////
    voidfractionNext_.correctBoundaryConditions();
	//for(int aa=0; aa< 1000; aa++)
	//{Info << voidfractionNext_[aa] << " voidfraction=" << endl; }	
    //用于修正或确保 voidfractionNext_ 数组中的值满足边界条件。这可能涉及到边界条件的处理、插值或其他纠正措施
    // reset counter of lost volume
    //if (verbose_) Pout << "Total particle volume neglected: " << tooMuch_<< endl; //输出体积分数被忽略的量，由于存在高于0.1的设置
    //tooMuch_ = 0.;
     if (verbose_)  	 
		 {
			 
			 
		 }
    // bring voidfraction from Eulerian Field to particle array
    //interpolationCellPoint<scalar> voidfractionInterpolator_(voidfractionNext_);
    //scalar voidfractionAtPos(0);

    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        label cellID = particleCloud_.cellIDs()[index][0];

        if(cellID >= 0)
        {
            voidfractions[index][0] = voidfractionNext_[cellID];
        }
        else
        {
            voidfractions[index][0] = -1.;
        }
    }
  
}

inline double CGVoidFraction::particleRadius(label index) const
{
    return particleCloud_.radius(index);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

