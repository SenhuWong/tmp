#include<cmath>
#define SR_EPS 0.000000000001
#define EQUAL(exp1,exp2) (fabs(exp1-exp2) < SR_EPS)
#define UNEQUAL(exp1,exp2) (fabs(exp1-exp2)>=SR_EPS) 
#define LESS(exp1,exp2) (exp1-exp2 <= -SR_EPS)
#define GREATER(exp1,exp2) (exp1-exp2 >=SR_EPS)
#define LEQUAL(exp1,exp2) (exp1-exp2 < SR_EPS)