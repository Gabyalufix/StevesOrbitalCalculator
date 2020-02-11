

lambert_problem_dTdx = function(DT,DDT,DDDT,x,T,lambda){
  L2 = lambda*lambda
  L3 = L2*lambda
  umx2 = 1 - x*x
  y = Math.sqrt(1 - L2*umx2)
  y2 = y*y
  y3 = y2*y
  return [
   1.0 / umx2 * (3.0 * T * x - 2.0 + 2.0 * L3 * x / y),
   1.0 / umx2 * (3.0 * T + 5.0 * x * DT + 2.0 * (1.0 - L2) * L3 / y3),
   1.0 / umx2 * (7.0 * x * DDT + 8.0 * DT - 6.0 * (1.0 - L2) * L2 * L3 * x / y3 / y2)
  ]
}


lambert_problem_x2tof2 = function(tof, x, N,lambda){
    a = 1.0 / (1.0 - x * x);
    if (a > 0) // ellipse
    {
        alfa = 2.0 * Math.acos(x);
        beta = 2.0 * Math.asin(Math.sqrt(lambda * lambda / a));
        if(lambda < 0.0) beta = -beta;
        tof = ((a * Math.sqrt(a) * ((alfa - Math.sin(alfa)) - (beta - Math.sin(beta)) + 2.0 * pi * N)) / 2.0);
		return tof;
    } else {
        alfa = 2.0 * Math.acosh(x);
        beta = 2.0 * Math.asinh(Math.sqrt(-lambda * lambda / a));
        if (lambda < 0.0) beta = -beta;
        tof = (-a * Math.sqrt(-a) * ((beta - Math.sinh(beta)) - (alfa - Math.sinh(alfa))) / 2.0);
		return tof;
    }
}

lambert_problem_hypergeometricF = function( z,  tol)
{
     Sj = 1.0;
     Cj = 1.0;
     err = 1.0;
     Cj1 = 0.0;
     Sj1 = 0.0;
     j = 0;
    while(err > tol) {
        Cj1 = Cj * (3.0 + j) * (1.0 + j) / (2.5 + j) * z / (j + 1);
        Sj1 = Sj + Cj1;
        err = Math.abs(Cj1);
        Sj = Sj1;
        Cj = Cj1;
        j = j + 1;
    }
    return Sj;
}

lambert_problem_x2tof = function(tof, x, N, lambda){
    var battin = 0.01;
    var lagrange = 0.2;
    var dist = Math.abs(x - 1);
    if (dist < lagrange && dist > battin) { // We use Lagrange tof expression
        return lambert_problem_x2tof2(tof, x, N, lambda);
    }
    var K = lambda * lambda;
    var E = x * x - 1.0;
    var rho = Math.abs(E);
    var z = Math.sqrt(1 + K * E);
    if (dist < battin) { // We use Battin series tof expression
        eta = z - lambda * x;
        S1 = 0.5 * (1.0 - lambda - x * eta);
        Q = lambert_problem_hypergeometricF(S1, 1e-11);
        Q = 4.0 / 3.0 * Q;
        tof = (eta * eta * eta * Q + 4.0 * lambda * eta) / 2.0 + N * pi / Math.pow(rho, 1.5);
        return tof;
    } else { // We use Lancaster tof expresion
        y = Math.sqrt(rho);
        g = x * z - lambda * E;
        d = 0.0;
        if (E < 0) {
            l = Math.acos(g);
            d = N * pi + l;
        } else {
            f = y * (z - lambda * x);
            d = Math.log(f + g);
        }
        tof = (x - lambda * z - d / y) / E;
        return tof;
    }
}





lambert_problem_householder = function( T, x0, N,  eps, iter_max, lambda){
     var it = 0;
     var err = 1.0;
     var xnew = 0.0;
     var tof = 0.0
	 var delta = 0.0
	 var DT = 0.0
	 var DDT = 0.0
	 var DDDT = 0.0;
    while((err > eps) && (it < iter_max)) {
        tof = lambert_problem_x2tof(tof, x0, N,lambda);
        dtdx = lambert_problem_dTdx(DT, DDT, DDDT, x0, tof,lambda);
		DT = dtdx[0]
		DDT = dtdx[1]
		DDDT = dtdx[2]
        delta = tof - T;
        DT2 = DT * DT;
        xnew = x0 - delta * (DT2 - delta * DDT / 2.0) / (DT * (DT2 - delta * DDT) + DDDT * delta * delta / 6.0);
        err = Math.abs(x0 - xnew);
        x0 = xnew;
        it++;
    }
    return [x0,it];
}





/** Constructs and solves a Lambert problem.
 *
 * \param[in] R1 first cartesian position
 * \param[in] R2 second cartesian position
 * \param[in] tof time of flight
 * \param[in] mu gravity parameter
 * \param[in] cw when 1 a retrograde orbit is assumed
 * \param[in] multi_revs maximum number of multirevolutions to compute
 */
lambert_problem_solve = function(r1, r2, tof, mu,
                                 cw, multi_revs){
    // 0 - Sanity checks
    if(tof <= 0) {
        console.log("Time of flight is negative!");
		return null;
    }
    if(mu <= 0) {
        console.log("Gravity parameter is zero or negative!");
		return null;
    }
    // 1 - Getting lambda and T
    m_c = Math.sqrt((r2[0] - r1[0]) * (r2[0] - r1[0]) + (r2[1] - r1[1]) * (r2[1] - r1[1])
               + (r2[2] - r1[2]) * (r2[2] - r1[2]));
     R1 = mag(r1);
     R2 = mag(r2);
    m_s = (m_c + R1 + R2) / 2.0;
	
	ir1 = vers(r1)
	ir2 = vers(r2)
	ih = cross(ir1,ir2)
	ih = vers(ih);
	
    if(ih[2] == 0){
        console.log("The angular momentum vector has no z component, impossible to define automatically clock or "+
                          "counterclockwise");
						  return null;
    }
    lambda2 = 1.0 - m_c / m_s;
    m_lambda = Math.sqrt(lambda2);
	lambda = m_lambda;

    if (ih[2] < 0.0) // Transfer angle is larger than 180 degrees as seen from abive the z axis
    {
        m_lambda = -m_lambda;
        it1 = cross(ir1, ih);
        it2 = cross(ir2, ih);
    } else {
        it1 = cross(ih, ir1);
        it2 = cross(ih, ir2);
    }
    it1 = vers(it1);
    it2 = vers(it2);

    if(cw){ // Retrograde motion
        m_lambda = -m_lambda;
        it1[0] = -it1[0];
        it1[1] = -it1[1];
        it1[2] = -it1[2];
        it2[0] = -it2[0];
        it2[1] = -it2[1];
        it2[2] = -it2[2];
    }
    lambda3 = m_lambda * lambda2;
    T = Math.sqrt(2.0 * mu / m_s / m_s / m_s) * tof;

    // 2 - We now have lambda, T and we will find all x
    // 2.1 - Let us first detect the maximum number of revolutions for which there exists a solution
    m_Nmax = Math.floor(T / pi);
    T00 = Math.acos(m_lambda) + m_lambda * Math.sqrt(1.0 - lambda2);
    T0 = (T00 + m_Nmax * pi);
    T1 = 2.0 / 3.0 * (1.0 - lambda3)
	DT = 0.0
	DDT = 0.0
	DDDT = 0.0;
    if(m_Nmax > 0) {
        if(T < T0) { // We use Halley iterations to find xM and TM
            it = 0;
            err = 1.0;
            T_min = T0;
            x_old = 0.0, x_new = 0.0;
            while(1) {
                dtdx = lambert_problem_dTdx(DT, DDT, DDDT, x_old, T_min,lambda);
				DT = dtdx[0]
				DDT = dtdx[1]
				DDDT = dtdx[2]
                if(DT != 0.0) {
                    x_new = x_old - DT * DDT / (DDT * DDT - DT * DDDT / 2.0);
                }
                err = Math.abs(x_old - x_new);
                if ((err < 1e-13) || (it > 12)) {
                    break;
                }
                tof = lambert_problem_x2tof(T_min, x_new, m_Nmax);
                x_old = x_new;
                it++;
            }
            if (T_min > T) {
                m_Nmax -= 1;
            }
        }
    }
    // We exit this if clause with Mmax being the maximum number of revolutions
    // for which there exists a solution. We crop it to m_multi_revs
    m_Nmax = Math.min(multi_revs, m_Nmax);

    // 2.2 We now allocate the memory for the output variables
    //m_v1.resize(m_Nmax * 2 + 1);
    //m_v2.resize(m_Nmax * 2 + 1);
    //m_iters.resize(m_Nmax * 2 + 1);
    //m_x.resize(m_Nmax * 2 + 1);
	 
	 m_x_size = m_Nmax * 2 + 1
     m_x = [];
	 m_v1 = []
	 m_v2 = []
	 m_iters = []
	 
    // 3 - We may now find all solutions in x,y
    // 3.1 0 rev solution
    // 3.1.1 initial guess
    if (T >= T00) {
        m_x[0] = -(T - T00) / (T - T00 + 4);
    } else if (T <= T1) {
        m_x[0] = T1 * (T1 - T) / (2.0 / 5.0 * (1 - lambda2 * lambda3) * T) + 1;
    } else {
        m_x[0] = Math.pow((T / T00), 0.69314718055994529 / Math.log(T1 / T00)) - 1.0;
    }
    // 3.1.2 Householder iterations
	hhold = lambert_problem_householder(T, m_x[0], 0, 1e-5, 15, lambda);
	m_x[0] = hhold[0];
    m_iters[0] = hhold[1]
    // 3.2 multi rev solutions
    var tmp;
    for( i = 1; i < m_Nmax + 1; i++) {
        // 3.2.1 left Householder iterations
        tmp = Math.pow((i * pi + pi) / (8.0 * T), 2.0 / 3.0);
        m_x[2 * i - 1] = (tmp - 1) / (tmp + 1);
		hhold = lambert_problem_householder(T, m_x[2 * i - 1], i, 1e-8, 15, lambda);
        m_iters[2 * i - 1] = hhold[1]
		m_x[2 * i - 1] = hhold[0];
        // 3.2.1 right Householder iterations
        tmp = Math.pow((8.0 * T) / (i * pi), 2.0 / 3.0);
        m_x[2 * i] = (tmp - 1) / (tmp + 1);
		hhold = lambert_problem_householder(T, m_x[2 * i], i, 1e-8, 15, lambda);
        m_iters[2 * i] = hhold[1];
		m_x[2 * i] = hhold[0]
    }

    // 4 - For each found x value we reconstruct the terminal velocities
     gamma = Math.sqrt(mu * m_s / 2.0);
     rho = (R1 - R2) / m_c;
     sigma =  Math.sqrt(1 - rho * rho);
    //double vr1, vt1, vr2, vt2, y;
    for( i = 0; i < m_x_size; i++) {
        y = Math.sqrt(1.0 - lambda2 + lambda2 * m_x[i] * m_x[i]);
        vr1 = gamma * ((m_lambda * y - m_x[i]) - rho * (m_lambda * y + m_x[i])) / R1;
        vr2 = -gamma * ((m_lambda * y - m_x[i]) + rho * (m_lambda * y + m_x[i])) / R2;
        vt = gamma * sigma * (y + m_lambda * m_x[i]);
        vt1 = vt / R1;
        vt2 = vt / R2;
		m_v1[i] = [];
		m_v2[i] = [];
        for( j = 0; j < 3; j++){
            m_v1[i][j] = vr1 * ir1[j] + vt1 * it1[j];
		}
        for( j = 0; j < 3; j++){
            m_v2[i][j] = vr2 * ir2[j] + vt2 * it2[j];
		}
    }
	var lambertSolution = {
		r1:r1,
		r2:r2,
		tof:tof,
		tofMonths : tof/onemonth,
		gamma:gamma,
		rho:rho,
		sigma:sigma,
		v1:m_v1,
		v2:m_v2,
		iters:m_iters
	}
	return lambertSolution;
}


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


/*
r1 = [1 * au,0,0];
r2 = [0,2*au,0];
tof = oneyear
mu = G * solmass
lambertSolve = lambert_problem_solve(r1, r2, tof,mu,0,0)
getOrbitalStatsFromRVM(r1,lambertSolve["v1"][0],mu)


lambertSolve

lambertSolve = lambert_problem_solve(r1, r2, tof/2,mu,0,0)
getOrbitalStatsFromRVM(r1,lambertSolve["v1"][0],mu)



v1 = lambertSolve["v1"][0]
v2 = lambertSolve["v2"][0]

cross( r1,lambertSolve["v1"] )

stats1 = getOrbitalStatsFromRVM(r1,v1,mu)
stats2 = getOrbitalStatsFromRVM(r2,v2,mu)

rv=r1
vv=v1
  hv = cross(rv,vv)
  nhat = cross([0,0,1],hv);
  v = mag(vv)
  r = mag(rv)
  h = mag(hv)

 vectorMult((vectorMult(r1,((mag(v)*mag(v) - mu/mag(r1)))) - vectorMult(v1,dot(r1,v1))) / mu


stats1

stats2
*/




