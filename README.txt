Input for the CG model

CGProps //Voidfraction
{
    alphaMin 0.1;
    scaleUpVol 1.0;
TopologyAccel true; Open Table-based method (only in Cartesian Grid!!!)
meshx X;//grid number in x-direction (only FOR TopologyAccel)
meshy Y;//grid number in y-direction (only FOR TopologyAccel)
meshz Z;//grid number in z-direction (only FOR TopologyAccel)
decx x;//processor number in z-direction (only FOR TopologyAccel)
decy y;//processor number in z-direction (only FOR TopologyAccel)
decz z;//processor number in z-direction (only FOR TopologyAccel)
}


CGlucy15DgidaspowProps //Force
{
    //verbose true;
    velFieldName "U";
    voidfractionFieldName "voidfraction";
    granVelFieldName "Us";
    //interpolation true;
    phi 1;
}