class DecayFit {
	public:	
		decay_parent(double *dim, double *par)
		{
			 /* function for fitting decay curves
			  a is amplitude of source,
			  b is tau of the source,
			  bkgrd is backround,
			  this first term is the exponetial decay
			 */
			 double a = par[0];
			 double b = par[1];
			 double x = dim[0];
			  
			 return a * TMath::Exp(-(b) * (x));
		}
		 double decay_daughter(double *dim,double *par) 
		 {  
			  /* function for fitting decay curves  
			   a is amplitude of source, 
			   b is tau of the source, 
			   c is tau of the daughter, 
			   bkgrd is backround, 
			    batemans equation. 
			   */ 
			   double a = par[0]; 
			   double b = par[1]; 
			   double c = par[2]; 
			   double d = par[3];  
			   double x = dim[0];
  
			   return d*((b*a)/(c-b))*(TMath::Exp(-(b)*(x))-TMath::Exp(-(c)*(x))); 
		 }
		 
