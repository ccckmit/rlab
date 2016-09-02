var R = require("../rlab");

// ODE (ordinary differential equation) : dopri() is an implementation of the Dormand-Prince-Runge-Kutta integrator with adaptive time-stepping
// https://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method
var ode = R.ODE(0,1,1,(t,y)=>y);
log('ode(0,1,1,(t,y)=>y)=', ode.strM());
log('at([0.3,0.7])=', ode.at(0.3,0.7));

// linear programming
log('LP:', R.solveLP([1,1],                   /* minimize [1,1]*x                */
                [[-1,0],[0,-1],[-1,-2]], /* matrix of inequalities          */
                [0,0,-3])                /* right-hand-side of inequalities */
    );

// Unconstrained minimization
var sqr = (x)=>x*x;
log('minimize(...) =', R.minimize((x)=>sqr(10*(x[1]-x[0]*x[0]))+sqr(1-x[0]),[-1.2,1]));