#include "Hungarian.h"

int main() {
	  /* an example cost matrix */
  int r[3*3] =  {20,90,10,60,30,40,90,90, 120};
  std::vector< std::vector<int> > m;
  int k = 0;

  m.resize(3);
  for (unsigned int i = 0; i < m.size(); i++)
  	m[i].resize(3);

  for (unsigned int i = 0; i < m.size(); i++)
  	for (unsigned int j = 0; j < m[i].size(); j++)
  		m[i][j] = r[k++];

  /* initialize the hungarian_problem using the cost matrix*/
  Hungarian hungarian(m , 3,3, HUNGARIAN_MODE_MAXIMIZE_UTIL) ;

  //fprintf(stderr, "assignement matrix has a now a size %d rows and %d columns.\n\n",  hungarian.ro,matrix_size);

  /* some output */
  fprintf(stderr, "cost-matrix:");
  hungarian.print_cost();

  /* solve the assignement problem */
  hungarian.solve();

  /* some output */
  fprintf(stderr, "assignment:");
  hungarian.print_assignment();

  std::map<std::pair<int, int>,int> foo;
  foo = hungarian.get_assignments();

  for (std::map<std::pair<int, int>,int>::iterator it = foo.begin(); it != foo.end(); ++it)
  	fprintf(stderr, "%d ",it->second);
  fprintf(stderr, "\n");


  return 0;
}