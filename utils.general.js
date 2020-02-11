
/**************************************************************************************
 * VECTOR ARITHMATIC:
 */

//dot product:
dot = function(x,y){
  return [ x[0]*y[0] + x[1]*y[1] + x[2]*y[2] ];
}

//cross-product
cross = function(rv,vv){
  return [ rv[1]*vv[2] - rv[2]*vv[1],
  rv[2]*vv[0] - rv[0]*vv[2],
  rv[0]*vv[1] - rv[1]*vv[0] ]
}
//magnitude:
mag = function(h){
  return Math.sqrt(h[0]*h[0] + h[1]*h[1] + h[2]*h[2])
}
//unit vector:
vers = function(v){
  //#gets the unit vector:
  vm = mag(v);
  return [v[0]/vm,v[1]/vm,v[2]/vm]
}
//multiply a vector times a constant:
vectorMult = function(v,k){
	return [ v[0]*k, v[1] * k, v[2]*k ]
}
//add two vectors:
vectorAdd = function(v1,v2){
	return [v1[0]+v2[0],v1[1]+v2[1],v1[2]+v2[2]];
}
