


getOrbitalStatsFromRVM = function(rv,vv,mu){
  hv = cross(rv,vv)
  nhat = cross([0,0,1],hv);
  v = mag(vv)
  r = mag(rv)
  h = mag(hv)
  ev = vectorMult( vectorAdd(vectorMult(rv,((v*v - mu/r))), vectorMult(vv,-dot(rv,vv))), 1/ mu)
  e = mag(ev)
  epsilon = (v*v / 2) - (mu / r)
  a = - 0.5 * mu / epsilon;
  p = a * (1 - e*e)
  i = Math.acos(hv[2]/h)
  var omega_longAscending
  var argp
  if(vv[2] == 0 && rv[2] == 0){
    omega_longAscending = 0;
    argp = Math.atan2( ev[1],ev[0] )
  } else {
    omega_longAscending = Math.acos(nhat[0]/mag(nhat))
    argp = Math.acos( dot(nhat,ev)/(mag(nhat)*e) )
  }
  if(nhat[1] < 0){
    omega_longAscending = 2*pi - omega_longAscending;
  }  
  if(ev[2] < 0){
    argp = 2*pi - argp
  }
  zz = dot(ev,rv)/( e*r )
  true_anomaly = Math.acos( Math.max(-1,Math.min(1,zz)) )
  if( dot(rv,vv) < 0){
    true_anomaly = 2*pi - true_anomaly
  }
  return {
    rv:rv,vv:vv,mu:mu,
    hv:hv,v:v,r:r,h:mag(hv),
    ev:ev,e:e,
    epsilon:  epsilon,
    a:a,p:p, 
    i:i,
	omega_longAscending:omega_longAscending,argp:argp,
	aInAU: a/au,
	rvInAU: vectorMult(rv,1/au),
    true_anomaly:true_anomaly
  }
}

addOrbitalStats = function(orbitstat,pstat){
	initialSize = 0 //Object.keys(stats).length
	
	while( initialSize != Object.keys(orbitstat).length + Object.keys(pstat).length){
		initialSize = Object.keys(orbitstat).length + Object.keys(pstat).length		
		
		if( orbitstat["b"] == null && orbitstat["a"] != null && orbitstat["e"] != null){
			orbitstat["b"] = orbitstat["a"] * (1-orbitstat["e"])
		}
		if( orbitstat["p"] == null && orbitstat["a"] != null && orbitstat["e"] != null){
			orbitstat["p"] = orbitstat["a"] * (1-orbitstat["e"]*orbitstat["e"])
		}
		if( orbitstat["epsilon"] == null && orbitstat["mu"] != null && orbitstat["a"] != null){
			orbitstat["epsilon"] = -orbitstat["mu"] / (2*orbitstat["a"])
		}
		if( orbitstat["epsilon"] == null && orbitstat["mu"] != null && pstat["v"] != null && pstat["r"] != null){
			orbitstat["epsilon"] = pstat["v"]*pstat["v"]/2 - orbitstat["mu"]/pstat["r"];
		}
		if( orbitstat["epsilon"] == null && orbitstat["mu"] != null && orbitstat["h"] != null && orbitstat["e"] != null){
			orbitstat["epsilon"] = -0.5 * ((orbitstat["mu"] * orbitstat["mu"]) / (orbitstat["h"]*orbitstat["h"]))*(1-orbitstat["e"]*orbitstat["e"])
		}
		
		//Positional stats:
		if( pstat["eccentric.anomaly"] == null && orbitstat["e"] != null && pstat["true.anomaly"] != null){
			if( orbitstat["e"] <= 1){
			   var costheta = Math.cos( pstat["true.anomaly"] )
			   pstat["eccentric.anomaly"] = Math.acos( ( orbitstat["e"] + costheta ) / ( 1 + orbitstat["e"] * costheta ) )
			   if( pstat["true.anomaly"] > pi ){
				   pstat["eccentric.anomaly"] = 2*pi - pstat["eccentric.anomaly"];
			   }
			} else {
			   var costheta = Math.cos( pstat["true.anomaly"] )
			   pstat["eccentric.anomaly"] = Math.acosh(( orbitstat["e"] + costheta ) / ( 1 + orbitstat["e"] * costheta ) )
			   if( pstat["true.anomaly"] > pi || pstat["true.anomaly"] < 0){
				   pstat["eccentric.anomaly"] = - pstat["eccentric.anomaly"];
			   }
			}
		}
		if( pstat["mean.anomaly"] == null && orbitstat["e"] != null && pstat["eccentric.anomaly"] != null){
			if( orbitstat["e"] <= 1){
			   pstat["mean.anomaly"] = pstat["eccentric.anomaly"] - orbitstat["e"] * Math.sin( pstat["eccentric.anomaly"] )
			} else {
               pstat["mean.anomaly"] = orbitstat["e"] * Math.sinh( pstat["eccentric.anomaly"] ) - pstat["eccentric.anomaly"]
			}
		}
		
	}
	return {orbitstat:orbitstat, pstat:pstat};
}

