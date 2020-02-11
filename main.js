



var canvas = document.getElementById("mainCanvas");
var ctx = canvas.getContext("2d");

var canvasHT = 400;
var canvasWT = 400;
ctx.lineWidth = 1;



drawXYcurve = function(xy,offset=200){
	x = xy["x"]
	y = xy["y"]
	ctx.moveTo( x[0]+offset,canvasHT-y[0]-offset);
	for(i =1; i<x.length;i++){
		ctx.lineTo(x[i]+offset, canvasHT-y[i]-offset);
	}
	ctx.stroke();
}
rescaleXY = function(xy,scalingFactor){
	rx = []
	ry = []
	x = xy["x"]
	y = xy["y"]
	for(i =0; i < x.length;i++){
		rx[i] = x[i] * scalingFactor;
		ry[i] = y[i] * scalingFactor;
	}
	return { x :rx, y:ry }
}


r1 = [1 * au,0,0];
r2 = [0,2*au,0];
tof = oneyear*1
mu = G * solmass
lambertSolve = lambert_problem_solve(r1, r2, tof,mu,0,0)


r1 = [1 * au,0,0];
r2 = [0,2*au,0];
tof = oneyear*0.50
mu = G * solmass
lambertSolve = lambert_problem_solve(r1, r2, tof,mu,0,0)
getOrbitalStatsFromRVM(r1,lambertSolve["v1"][0],mu)


lambertSolve = lambert_problem_solve(r1, r2, tof*0.1,mu,0,0)
getOrbitalStatsFromRVM(r1,lambertSolve["v1"][0],mu)

lambertSolve = lambert_problem_solve(r1, r2, tof*1.5,mu,0,0)
getOrbitalStatsFromRVM(r1,lambertSolve["v1"][0],mu)

//ctx.moveTo(0, 0);
//ctx.lineTo(200, 100);
//ctx.stroke();
/*
r1 = [1 * au,0,0];
r2 = [0,2*au,0];
tof = oneyear
mu = G * solmass
lambertSolve = lambert_problem_solve(r1, r2, tof,mu,0,0)
lambertSolve

v1 = lambertSolve["v1"][0]
v2 = lambertSolve["v2"][0]

cross( r1,lambertSolve["v1"] )

stats1 = getOrbitalStatsFromRVM(r1,v1,mu)
stats2 = getOrbitalStatsFromRVM(r2,v2,mu)

aInAU = stats1["aInAU"]
e = stats1["e"]

ept = getEllipsePt(a,e,0,100)
ept = rescaleXY(ept, 10/au)

x = ept["x"]
y = ept["y"]


ctx.lineTo(x[50], canvasHT-y[50])
ctx.stroke();

drawXYcurve(ept)

//getEllipsePt
ctx.moveTo( x[0]+offset,canvasHT-y[0]-offset);
ctx.lineTo(x[50]+offset, canvasHT-y[50]-offset)
ctx.stroke();


[x[0]+offset, canvasHT-y[0]-offset]
[x[50]+offset, canvasHT-y[50]-offset]



ctx.moveTo(0, 0);
ctx.lineTo(200, 100);
ctx.lineTo(100, 200);
ctx.stroke();*/