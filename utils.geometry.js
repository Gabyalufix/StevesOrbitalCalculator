

getEllipsePt = function(a,e,offset=0,pointct=10000){
 p = a * (1 - e^2);
 x = []
 y = []
 r = []
 tt = []
 t  = []
 for( i = 0; i < pointct; i++){
	 tt[i] = (i / (pointct-1)) * 2 * pi;
	 r[i] = p / ( 1 + e*Math.cos(tt[i]));
	 x[i] = r[i] * Math.cos(tt[i]+offset);
	 y[i] = r[i] * Math.sin(tt[i]+offset);
 }
 return {
	x:x,y:y,r:r,
	theta : tt,
    t:t,a:a,e:e,p:p,offset:offset,pointct:pointct
 };
}
