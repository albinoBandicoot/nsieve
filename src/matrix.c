#include "matrix.h"

void solve_matrix (nsieve_t *ns){
}

/* Factor deduction and dealing with multiple polynomials, etc.
 *
 * First, note that we are trying to produce congruences of the form X^2 ~= Y^2 (mod N). gcd (X-Y, N) will (about half
 * of the time) yield a nontrivial factor of N. Each polynomial Q(x) is constructed so that 
 *
 * 	a Q(x) = (ax + b)^2 - N   ~= H^2 (mod N), for H = (ax + b) (mod n) 
 *
 * The sieve will find many values of x such that Q(x) = y is smooth over the factor base. Note that a*Q(x) may or may
 * not actually factor over the factor base (this depends on whether the g_i are in the factor base or not, and we cannot
 * assume that they will be), but we really need the LHS to be a square mod N, so we must keep this factor a. This creates
 * a puzzle, until we realize that we can select one relation, say   a * Q(x_0) = y_0   and multiply all of the other 
 * relations that share the same 'a' value by it. We then have a list of relations like this:
 *
 * 	a^2 Q(x_1) Q(x_0) = y_1 * y_0
 * 	...
 * 	a^2 Q(x_n) Q(x_0) = y_n * y_0
 *
 * Now all of the LHS's are squares mod N, and the RHS's are still smooth over the factor base. Define
 *
 * 	H_q,i = [(ax_i + b) * (ax_0 + b)] % N	 for the a,b values of the polynomial q. 
 *
 * Then we may rewrite the relations above like this:
 *
 * 	(H_Q,i)^2 ~= y_i * y_0   (mod N).
 *
 * The linear algebra step will then select a subset of these relations so that the RHS is a square (in the integers). 
 * Note that the relations may come from different polynomials. We get
 *
 * 	(H_p,i)^2 * (H_q,j)^2 * ... ~= y_p,i * y_q,j ...  (mod N)
 * 	   "           "         "  ~= v^2		  (mod N)
 *
 * This is our desired congruence of squres. 
*/

void deduce_factors (nsieve_t *ns /* and history matrix */){
}
