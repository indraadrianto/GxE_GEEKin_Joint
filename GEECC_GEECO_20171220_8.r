#updates on 12/20/2017
#add steps to allow coVariates in GEE_Cpp_CO function
#add steps to allow missing genotypes in GEE_Cpp_CO function

#add steps to allow missing genotypes with GEECC_GEECO_20170330_6.R.R file




#GEE model with kinship as the working correlation matrix
#y: dependent variable (a binary variable eg. disease: yes or no)
#XC: independent variables such as environmental factors
#fid: family ids
#K: n*n kinship matrix (its the pairwise correlation between subjects)

#------------fit a QLS ---------------#



########## Test R-C++ integration ############
library(Rcpp)
library(RcppArmadillo)


cppFunction('void  GEE_Cpp_CC(const char* OutFile, const char* ErrorLog, arma::mat & BI, arma::vec & Y0 , arma::vec & E, arma::mat & SNP, arma::mat & CoVar,arma::vec & fam0, arma::mat & K0, arma::mat & Kinv0, int threads=1){

	//run the model y ~ b0+e to estimate dnew once
	
		//std::clock_t start1 = clock();

		
		double tol = 0.000001;
	
	//if( threads > 0 ){
	//omp_set_num_threads( threads );
	//}
	//std::cout<<"cores:"<<omp_get_max_threads()<<std::endl;

//arma::mat outcome= arma::zeros<arma::mat>(2*BI.n_cols,3);
std::ofstream output;
output.open(OutFile);

std::ofstream errorouput;
errorouput.open(ErrorLog);

output <<"snp.id"<<"\t"<< "coef"<<"\t"<<"sd"<<"\t"<<"p-value"<<"\t"<<"Z-value"<<"\\n";

int nvars = 4 + CoVar.n_cols;

//#pragma omp parallel for schedule(dynamic)
for(int sp = 0; sp<BI.n_cols; sp++){
	//std::cout<<"Start processing SNP:"<<sp<<std::endl;
	
	//std::clock_t start2 = clock();
	
	arma::vec snps = SNP.col(sp);
	int er = 0;
	double d = 0.5;
	double dnew = d;
	double diff1 = 1;
	double diff2 = diff1;
	double tol = 0.00001;
	int niter = 0;
	int maxiter = 100;
	//if(!is.finite(snps)){
		arma::uvec finiteids = arma::find_finite(snps);
		int ny = finiteids.n_elem;
	//std::cout<<"number of nonmissing is "<<ny<<"\\n"<<std::endl;

		arma::mat XC = arma::ones<arma::mat>(ny,nvars);
		XC.col(1) = snps(finiteids);
		XC.col(2) = E(finiteids);
		XC.col(3) = XC.col(1) % XC.col(2);
		XC.cols(4,nvars-1) = CoVar.rows(finiteids);
		
		arma:: vec Y = Y0(finiteids);
		arma::vec bini = BI.col(sp);
		arma::vec b = bini;
		arma::mat tempdat = arma::zeros<arma::mat>(ny,ny);
		arma::mat D = arma::zeros<arma::mat>(ny,nvars);
		arma::mat Vinv = tempdat;
		arma::mat Dv = arma::zeros<arma::mat>(nvars,ny);
		arma::mat DVDinv = arma::zeros<arma::mat>(nvars,nvars);
		arma::vec pipi(ny);
		arma::mat A = arma::zeros<arma::mat>(ny,ny);
		arma::mat K = K0(finiteids,finiteids);
		arma::mat Kinv = Kinv0(finiteids,finiteids);
		arma::vec fam = fam0(finiteids);
		
  		std::map<int, int> mapcounts;

 		int n = fam.size();
 		for (int i = 0; i < n; i++) {
 		   mapcounts[fam[i]]++;
 		}
		
		arma::vec ni(mapcounts.size());
		int i=0;

		std::map<int, int>::iterator it = mapcounts.begin();
    		while(it != mapcounts.end())
    		{
		    ni(i)=it->second;
    		    it++;
		    i++;
   		 }

	//std::cout<<"number of ni is "<<ni.size()<<"\\n"<<std::endl;

	


	//int ny = Y.n_elem;
	//arma::mat XC = arma::ones<arma::mat>(ny,nvars);
	//XC.ones;
//the order of the colum is important and has to be consistant with initial values
	//XC.col(1) = SNP.col(sp);
	//XC.col(2) = E;
	//XC.col(3) = XC.col(1) % XC.col(2);
	//XC.cols(4,nvars-1) = CoVar;

	//arma::vec bini = BI.col(sp);
	//arma::vec b = bini;
	//arma::mat tempdat = arma::zeros<arma::mat>(ny,ny);
	//arma::mat D = arma::zeros<arma::mat>(ny,nvars);
	//arma::mat Vinv = tempdat;
	//arma::mat Dv = arma::zeros<arma::mat>(nvars,ny);
	//arma::mat DVDinv = arma::zeros<arma::mat>(nvars,nvars);
	//arma::vec pipi(ny);
	//arma::mat A = arma::zeros<arma::mat>(ny,ny);
	//std::cout<<"sucessfull here"<<std::endl;

	//std::clock_t end2 = clock();
	//double time2 = (double)(end2-start2)/CLOCKS_PER_SEC;
	//std::cout<<"time before while loop is "<<time2<<"\\n"<<std::endl;


	//while(diff1 > tol && diff2 > tol && niter < maxiter){
	while((diff1 > tol || diff2 > tol) &&  niter < maxiter){
	
	//std::clock_t start20 = clock();

	arma::vec ebx = exp(XC * bini); //may need a loop here to do every elements of ebx
	arma::vec pi = ebx/(1.0+ebx);
	

	//std::clock_t start200 = clock();
	//double time211 = (double)(start200-start20)/CLOCKS_PER_SEC;
	//std::cout<<"time calculating block1_1 is "<<time211<<"\\n"<<std::endl;
	
	//arma::mat A = arma::diagmat(pi % (1-pi)); // this is element wise production
	//arma::mat B = arma::diagmat(pi);
	
	//A.diag() = pipi;

	for(int ii=0;ii<ny;ii++){
		pipi(ii) = pi(ii)*(1-pi(ii));
		A(ii,ii) = pipi(ii);
	}

	//std::clock_t start201 = clock();
	//double time212 = (double)(start201-start200)/CLOCKS_PER_SEC;
	//std::cout<<"time calculating block1_2 is "<<time212*1000<<"\\n"<<std::endl;

	//D = A * XC;
	//use the fact that A is a diagnal matrix to simplify the calculation
	for(int ii=0;ii<nvars;ii++){
		for(int jj=0;jj<ny;jj++){
			D(jj,ii) = A(jj,jj) * XC(jj,ii);
		}
	}
	
	//std::clock_t start202 = clock();
	//double time213 = (double)(start202-start201)/CLOCKS_PER_SEC;
	//std::cout<<"time calculating block1_3 is "<<time213<<"\\n"<<std::endl;

	
	NumericVector piip = as<NumericVector>(wrap(pipi));
	NumericVector piipsq = sqrt(1/piip);
	arma::vec Asqiv = as<arma::vec>(piipsq);

	//std::clock_t start203 = clock();
	//double time214 = (double)(start203-start202)/CLOCKS_PER_SEC;
	//std::cout<<"time calculating block1_4 is "<<time214<<"\\n"<<std::endl;


	//std::cout<<"sucessfull here 2"<<std::endl;

	//Vinv = ( Asqiv * arma::trans(Asqiv) ) % Kinv * 1/d;
	//for(int ii = 0; ii<Vinv.n_cols; ii++){
	//	Vinv.col(ii) = Asqiv(ii) * Asqiv % Kinv.col(ii) * 1/d;
	//}
	int st = 0;
	int ed = -1;
	int l = 0;
	arma::rowvec TAsqiv = arma::trans(Asqiv);
		

	//std::cout<<"TAsqinv is sucessful"<<std::endl;
	
	//std::clock_t start21 = clock();
	
	//double time21 = (double)(start21-start20)/CLOCKS_PER_SEC;
	//std::cout<<"time calculating block1 is "<<time21<<"\\n"<<std::endl;

	arma::mat Dtr = arma::trans(D);
	
	for(int ii=0;ii<ni.n_elem;ii++){
		l = ni(ii);
		ed = ed+l;
		if (l == 1) {
			Vinv(ed,ed) = Asqiv(ed) * TAsqiv(ed);
			Dv.col(ed) = Dtr.col(ed) * Vinv(ed,ed);
			  }
		else if (l > 1) {
			arma::mat R = K(arma::span(st,ed),arma::span(st,ed)) * d;
			R.diag().ones();
			arma::mat Rinv = inv(R);
			Vinv(arma::span(st,ed),arma::span(st,ed)) = Asqiv(arma::span(st,ed)) * TAsqiv(arma::span(st,ed)) % Rinv;
			Dv.cols(st,ed) = Dtr.cols(st,ed)*Vinv(arma::span(st,ed),arma::span(st,ed));
			}
		st = st + l;
	}

	//std::cout<<"Vinv is sucessful"<<std::endl;
	
	//Vinv.diag().ones();
	
	//std::clock_t start22 = clock();
	
	//double time22 = (double)(start22-start21)/CLOCKS_PER_SEC;
	//std::cout<<"time calculating block2 is "<<time22<<"\\n"<<std::endl;

	//Dv = arma::trans(D) * Vinv;
	//Vinv is a block matrix and use this to reduce computation burden
	
	
	//std::clock_t start231 = clock();
	//double time231 = (double)(start231-start22)/CLOCKS_PER_SEC;
	//std::cout<<"time calculating block3_1 is "<<time231<<"\\n"<<std::endl;
	
	arma::mat U = Dv * (Y-pi);
	arma::mat I = Dv * D;

	//std::clock_t start23 = clock();
	//double time23 = (double)(start23-start22)/CLOCKS_PER_SEC;
	//std::cout<<"time calculating block3 is "<<time23<<"\\n"<<std::endl;


	//std::cout<<"before inv(I) is sucessful"<<"\\n"<<std::endl;
	//std::cout<<I(0,0)<<","<<I(0,1)<<"\\n"<<std::endl;
	//std::cout<<I(1,0)<<","<<I(1,1)<<"\\n"<<std::endl;

	try
	{
		DVDinv = inv(I);
	}
	catch (...)
  	{
		er = 1;
		break;
 	}

	//std::cout<<"after inv(I) is sucessful"<<"\\n"<<std::endl;

	//std::clock_t start24 = clock();
	//double time24 = (double)(start24-start23)/CLOCKS_PER_SEC;
	//std::cout<<"time calculating block4 is "<<time24<<"\\n"<<std::endl;


	//update b
	b = bini + DVDinv * U;
	//diff1 = max(abs(1-b/bini)); //may need to break this
	arma::vec bt = abs(1-b/bini);
	diff1 = bt.max();
	bini = b;

	NumericVector XCb = as<NumericVector>(wrap(XC * b));
	NumericVector pinew = plogis(XCb);
	arma::vec Ypi = Y - as<arma::vec>(pinew);
	NumericVector pinewsq = sqrt(pinew*(1-pinew));
	arma::vec pisq = as<arma::vec>(pinewsq);
	//arma::vec rhat = Ypi/pisq;
	arma::vec rhat = Ypi;
	//std::cout<<"number of ni is "<<ni.size()<<"\\n"<<std::endl;


	st = 0;
	ed = -1;
	double mynu = 0;
	double myde = 0;
	l = 0;
	
	for(int i=0;i<ni.n_elem;i++){
		l = ni(i);
		ed = ed+l;
		if (l == 1) {
		    tempdat(ed,ed) = rhat(ed) * rhat(ed) ;
			  }
		else if (l > 1) {
		  	tempdat(arma::span(st,ed),arma::span(st,ed)) = rhat(arma::span(st,ed))* arma::trans(rhat(arma::span(st,ed)));
		        mynu = mynu + arma::accu(arma::trimatu(tempdat(arma::span(st,ed),arma::span(st,ed)))) - arma::trace(tempdat(arma::span(st,ed),arma::span(st,ed)));
		        myde = myde + arma::accu(arma::trimatu(K(arma::span(st,ed),arma::span(st,ed)))) - arma::trace(K(arma::span(st,ed),arma::span(st,ed)));

			  }
 		st = st + l;
	}
	dnew = mynu/myde;
	diff2 = std::abs(dnew -d);
	d = dnew;
	niter = niter+1;	


	} //last bracket for while loop;
	

	//std::clock_t end3 = clock();
	//double time3 = (double)(end3-end2)/CLOCKS_PER_SEC;
	//std::cout<<"niter is "<<niter<<"\\n"<<std::endl;
	//std::cout<<"time excuting while loop is "<<time3<<"\\n"<<std::endl;
	
	//std::cout<<"b is"<<b(0)<<" and "<<b(1)<<"\\n"<<std::endl;

	if(er==1) {
		errorouput<< "There is an error with SNP : " << sp << "\\n";
		continue;
	}
	
	
	//std::clock_t end31 = clock();
	//double time31 = (double)(end31-end3)/CLOCKS_PER_SEC;
	//std::cout<<"time excuting while loop 1_1 is "<<time31<<"\\n"<<std::endl;
	
	//use the properties of block matrix of tempdat and Vinv to increase calculation speed
	//arma::mat Vb =  DVDinv * Dv * tempdat * Vinv * D * DVDinv; // robust standard error of b
	//Vbtemp = Dv * tempdat * Vinv
	arma::mat Vbtemp = arma::zeros<arma::mat>(nvars,ny);
	
	int st = 0;
	int ed = -1;
	int l = 0;
	for(int ii=0;ii<ni.n_elem;ii++){
			l = ni(ii);
			ed = ed+l;
			if (l == 1) {
				Vbtemp.col(ed) = Dv.col(ed)*tempdat(ed,ed)*Vinv(ed,ed);
				  }
			else if (l > 1) {
				Vbtemp.cols(st,ed)=Dv.cols(st,ed)*tempdat(arma::span(st,ed),arma::span(st,ed))*Vinv(arma::span(st,ed),arma::span(st,ed));
				}
			st = st + l;
	}
	
	arma::mat Vb = DVDinv * Vbtemp * D * DVDinv;
	
	//std::clock_t end32 = clock();
	//double time32 = (double)(end32-end31)/CLOCKS_PER_SEC;
	//std::cout<<"time excuting while loop 1_2 is "<<time32<<"\\n"<<std::endl;

	NumericVector resb = as<NumericVector>(wrap(bini));
	arma::vec Vbdiag = Vb.diag();
	//std::cout<<"Vbdiag is "<<Vbdiag(0)<<Vbdiag(1)<<"\\n"<<std::endl;

	NumericVector resseb = sqrt(as<NumericVector>(wrap(Vbdiag)));
	NumericVector  resz = resb/resseb;
	NumericVector  p = pnorm(resz);
 	NumericVector resp = 2*pmin(p,1-p);

 	//std::stringstream buf;
	//for(int aa=0; aa<2; aa++){
	//buf <<sp<<"\t"<< bini(aa)<<"\t"<<resseb(aa)<<"\t"<<resp(aa)<<"\\n";
	//}

	//#pragma omp critical
	//output << buf.rdbuf();
 	//{
	for(int aa=0; aa<nvars; aa++){
		output <<sp<<"\t"<< resb(aa)<<"\t"<<resseb(aa)<<"\t"<<resp(aa)<<"\t"<<resz(aa)<<"\\n";
		}
	//}
	//std::cout<<"End Processing SNP:"<<sp<<std::endl;
	
	//std::clock_t end4 = clock();
	//double time4 = (double)(end4-end32)/CLOCKS_PER_SEC;
	//std::cout<<"time after excuting while loop is "<<time4<<"\\n"<<std::endl;

}

//return outcome;

output.close();

}',depends="RcppArmadillo", plugins = "openmp")











cppFunction('void  GEE_Cpp_CO(const char* OutFile, const char* ErrorLog, arma::mat & BI, arma::vec & Y0 , arma::mat & SNP, arma::mat & CoVar, arma::vec & fam0, arma::mat & K0, arma::mat & Kinv0, int threads=1){


std::ofstream output;
output.open(OutFile);

std::ofstream errorouput;
errorouput.open(ErrorLog);

output <<"snp.id"<<"\t"<< "coef"<<"\t"<<"sd"<<"\t"<<"p-value"<<"\t"<<"Z-value"<<"\\n";

int nvars = 2 + CoVar.n_cols;

//#pragma omp parallel for schedule(dynamic)
for(int sp = 0; sp<BI.n_cols; sp++){
	//std::cout<<"Start processing SNP:"<<sp<<std::endl;
	
	//std::clock_t start2 = clock();
	

	//int ny = Y.n_elem;
	arma::vec snps = SNP.col(sp);
	int er = 0;
	double d = 0.5;
	double dnew = d;
	double diff1 = 1;
	double diff2 = diff1;
	double tol = 0.00001;
	int niter = 0;
	int maxiter = 100;

	arma::uvec finiteids = arma::find_finite(snps);
	int ny = finiteids.n_elem;
	arma::mat XC = arma::ones<arma::mat>(ny,nvars);
	//XC.ones;
	//the order of the colum is important and has to be consistant with initial values
	XC.col(1) = snps(finiteids);
	XC.cols(2,nvars-1) = CoVar.rows(finiteids);
	arma:: vec Y = Y0(finiteids);

	arma::vec bini = BI.col(sp);
	arma::vec b = bini;
	arma::mat tempdat = arma::zeros<arma::mat>(ny,ny);
	arma::mat D = arma::zeros<arma::mat>(ny,nvars);
	arma::mat Vinv = tempdat;
	arma::mat Dv = arma::zeros<arma::mat>(nvars,ny);
	arma::mat DVDinv = arma::zeros<arma::mat>(nvars,nvars);
	arma::vec pipi(ny);
	arma::mat A = arma::zeros<arma::mat>(ny,ny);
	arma::mat K = K0(finiteids,finiteids);
	arma::mat Kinv = Kinv0(finiteids,finiteids);
	arma::vec fam = fam0(finiteids);

	std::map<int, int> mapcounts;

 	int n = fam.size();
 	for (int i = 0; i < n; i++) {
 	   mapcounts[fam[i]]++;
 	}
		
	arma::vec ni(mapcounts.size());
	int i=0;

	std::map<int, int>::iterator it = mapcounts.begin();
    	while(it != mapcounts.end())
    	{
	    ni(i)=it->second;
    	    it++;
	    i++;
   	 }

	//std::cout<<"sucessfull here"<<std::endl;

	//std::clock_t end2 = clock();
	//double time2 = (double)(end2-start2)/CLOCKS_PER_SEC;
	//std::cout<<"time before while loop is "<<time2<<"\\n"<<std::endl;


	//while(diff1 > tol && diff2 > tol && niter < maxiter){
	while((diff1 > tol || diff2 > tol) &&  niter < maxiter){
	
	//std::clock_t start20 = clock();

	arma::vec ebx = exp(XC * bini); //may need a loop here to do every elements of ebx
	arma::vec pi = ebx/(1.0+ebx);
	
	
	//std::clock_t start200 = clock();
	//double time211 = (double)(start200-start20)/CLOCKS_PER_SEC;
	//std::cout<<"time calculating block1_1 is "<<time211<<"\\n"<<std::endl;
	
	//arma::mat A = arma::diagmat(pi % (1-pi)); // this is element wise production
	//arma::mat B = arma::diagmat(pi);
	
	//A.diag() = pipi;

	for(int ii=0;ii<ny;ii++){
		pipi(ii) = pi(ii)*(1-pi(ii));
		A(ii,ii) = pipi(ii);
	}

	//std::clock_t start201 = clock();
	//double time212 = (double)(start201-start200)/CLOCKS_PER_SEC;
	//std::cout<<"time calculating block1_2 is "<<time212*1000<<"\\n"<<std::endl;

	//D = A * XC;
	//use the fact that A is a diagnal matrix to simplify the calculation
	for(int ii=0;ii<nvars;ii++){
		for(int jj=0;jj<ny;jj++){
			D(jj,ii) = A(jj,jj) * XC(jj,ii);
		}
	}
	
	//std::clock_t start202 = clock();
	//double time213 = (double)(start202-start201)/CLOCKS_PER_SEC;
	//std::cout<<"time calculating block1_3 is "<<time213<<"\\n"<<std::endl;

	
	NumericVector piip = as<NumericVector>(wrap(pipi));
	NumericVector piipsq = sqrt(1/piip);
	arma::vec Asqiv = as<arma::vec>(piipsq);

	//std::clock_t start203 = clock();
	//double time214 = (double)(start203-start202)/CLOCKS_PER_SEC;
	//std::cout<<"time calculating block1_4 is "<<time214<<"\\n"<<std::endl;


	//std::cout<<"sucessfull here 2"<<std::endl;

	//Vinv = ( Asqiv * arma::trans(Asqiv) ) % Kinv * 1/d;
	//for(int ii = 0; ii<Vinv.n_cols; ii++){
	//	Vinv.col(ii) = Asqiv(ii) * Asqiv % Kinv.col(ii) * 1/d;
	//}
	int st = 0;
	int ed = -1;
	int l = 0;
	arma::rowvec TAsqiv = arma::trans(Asqiv);

	//std::cout<<"TAsqinv is sucessful"<<std::endl;
	
	//std::clock_t start21 = clock();
	
	//double time21 = (double)(start21-start20)/CLOCKS_PER_SEC;
	//std::cout<<"time calculating block1 is "<<time21<<"\\n"<<std::endl;

	arma::mat Dtr = arma::trans(D);
	
	for(int ii=0;ii<ni.n_elem;ii++){
		l = ni(ii);
		ed = ed+l;
		if (l == 1) {
			Vinv(ed,ed) = Asqiv(ed) * TAsqiv(ed);
			Dv.col(ed) = Dtr.col(ed) * Vinv(ed,ed);
			  }
		else if (l > 1) {
			arma::mat R = K(arma::span(st,ed),arma::span(st,ed)) * d;
			R.diag().ones();
			arma::mat Rinv = inv(R);
			Vinv(arma::span(st,ed),arma::span(st,ed)) = Asqiv(arma::span(st,ed)) * TAsqiv(arma::span(st,ed)) % Rinv;
			Dv.cols(st,ed) = Dtr.cols(st,ed)*Vinv(arma::span(st,ed),arma::span(st,ed));
			}
		st = st + l;
	}

	//std::cout<<"Vinv is sucessful"<<std::endl;
	
	//Vinv.diag().ones();
	
	//std::clock_t start22 = clock();
	
	//double time22 = (double)(start22-start21)/CLOCKS_PER_SEC;
	//std::cout<<"time calculating block2 is "<<time22<<"\\n"<<std::endl;

	//Dv = arma::trans(D) * Vinv;
	//Vinv is a block matrix and use this to reduce computation burden
	
	
	//std::clock_t start231 = clock();
	//double time231 = (double)(start231-start22)/CLOCKS_PER_SEC;
	//std::cout<<"time calculating block3_1 is "<<time231<<"\\n"<<std::endl;
	
	arma::mat U = Dv * (Y-pi);
	arma::mat I = Dv * D;

	//std::clock_t start23 = clock();
	//double time23 = (double)(start23-start22)/CLOCKS_PER_SEC;
	//std::cout<<"time calculating block3 is "<<time23<<"\\n"<<std::endl;


	//std::cout<<"before inv(I) is sucessful"<<"\\n"<<std::endl;
	//std::cout<<I(0,0)<<","<<I(0,1)<<"\\n"<<std::endl;
	//std::cout<<I(1,0)<<","<<I(1,1)<<"\\n"<<std::endl;

	try
	{
		DVDinv = inv(I);
	}
	catch (...)
  	{
		er = 1;
		break;
 	}

	//std::cout<<"after inv(I) is sucessful"<<"\\n"<<std::endl;

	//std::clock_t start24 = clock();
	//double time24 = (double)(start24-start23)/CLOCKS_PER_SEC;
	//std::cout<<"time calculating block4 is "<<time24<<"\\n"<<std::endl;


	//update b
	b = bini + DVDinv * U;
	//diff1 = max(abs(1-b/bini)); //may need to break this
	arma::vec bt = abs(1-b/bini);
	diff1 = bt.max();
	bini = b;

	NumericVector XCb = as<NumericVector>(wrap(XC * b));
	NumericVector pinew = plogis(XCb);
	arma::vec Ypi = Y - as<arma::vec>(pinew);
	NumericVector pinewsq = sqrt(pinew*(1-pinew));
	arma::vec pisq = as<arma::vec>(pinewsq);
	//arma::vec rhat = Ypi/pisq;
	arma::vec rhat = Ypi;


	st = 0;
	ed = -1;
	double mynu = 0;
	double myde = 0;
	l = 0;
	
	for(int i=0;i<ni.n_elem;i++){
		l = ni(i);
		ed = ed+l;
		if (l == 1) {
		    tempdat(ed,ed) = rhat(ed) * rhat(ed) ;
			  }
		else if (l > 1) {
		  	tempdat(arma::span(st,ed),arma::span(st,ed)) = rhat(arma::span(st,ed))* arma::trans(rhat(arma::span(st,ed)));
		        mynu = mynu + arma::accu(arma::trimatu(tempdat(arma::span(st,ed),arma::span(st,ed)))) - arma::trace(tempdat(arma::span(st,ed),arma::span(st,ed)));
		        myde = myde + arma::accu(arma::trimatu(K(arma::span(st,ed),arma::span(st,ed)))) - arma::trace(K(arma::span(st,ed),arma::span(st,ed)));

			  }
 		st = st + l;
	}
	dnew = mynu/myde;
	diff2 = std::abs(dnew -d);
	d = dnew;
	niter = niter+1;	


	} //last bracket for while loop;
	
	//std::clock_t end3 = clock();
	//double time3 = (double)(end3-end2)/CLOCKS_PER_SEC;
	//std::cout<<"niter is "<<niter<<"\\n"<<std::endl;
	//std::cout<<"time excuting while loop is "<<time3<<"\\n"<<std::endl;
	
	//std::cout<<"b is"<<b(0)<<" and "<<b(1)<<"\\n"<<std::endl;

	if(er==1) {
		errorouput<< "There is an error with SNP : " << sp << "\\n";
		continue;
	}
	
	
	//std::clock_t end31 = clock();
	//double time31 = (double)(end31-end3)/CLOCKS_PER_SEC;
	//std::cout<<"time excuting while loop 1_1 is "<<time31<<"\\n"<<std::endl;
	
	//use the properties of block matrix of tempdat and Vinv to increase calculation speed
	//arma::mat Vb =  DVDinv * Dv * tempdat * Vinv * D * DVDinv; // robust standard error of b
	//Vbtemp = Dv * tempdat * Vinv
	arma::mat Vbtemp = arma::zeros<arma::mat>(nvars,ny);
	
	int st = 0;
	int ed = -1;
	int l = 0;
	for(int ii=0;ii<ni.n_elem;ii++){
			l = ni(ii);
			ed = ed+l;
			if (l == 1) {
				Vbtemp.col(ed) = Dv.col(ed)*tempdat(ed,ed)*Vinv(ed,ed);
				  }
			else if (l > 1) {
				Vbtemp.cols(st,ed)=Dv.cols(st,ed)*tempdat(arma::span(st,ed),arma::span(st,ed))*Vinv(arma::span(st,ed),arma::span(st,ed));
				}
			st = st + l;
	}
	
	arma::mat Vb = DVDinv * Vbtemp * D * DVDinv;
	
	//std::clock_t end32 = clock();
	//double time32 = (double)(end32-end31)/CLOCKS_PER_SEC;
	//std::cout<<"time excuting while loop 1_2 is "<<time32<<"\\n"<<std::endl;

	NumericVector resb = as<NumericVector>(wrap(bini));
	arma::vec Vbdiag = Vb.diag();
	//std::cout<<"Vbdiag is "<<Vbdiag(0)<<Vbdiag(1)<<"\\n"<<std::endl;

	NumericVector resseb = sqrt(as<NumericVector>(wrap(Vbdiag)));
	NumericVector  resz = resb/resseb;
	NumericVector  p = pnorm(resz);
 	NumericVector resp = 2*pmin(p,1-p);

 	//std::stringstream buf;
	//for(int aa=0; aa<2; aa++){
	//buf <<sp<<"\t"<< bini(aa)<<"\t"<<resseb(aa)<<"\t"<<resp(aa)<<"\\n";
	//}

	//#pragma omp critical
	//output << buf.rdbuf();
 	//{
	for(int aa=0; aa<nvars; aa++){
		output <<sp<<"\t"<< resb(aa)<<"\t"<<resseb(aa)<<"\t"<<resp(aa)<<"\t"<<resz(aa)<<"\\n";
		}
	//}
	//std::cout<<"End Processing SNP:"<<sp<<std::endl;
	
	//std::clock_t end4 = clock();
	//double time4 = (double)(end4-end32)/CLOCKS_PER_SEC;
	//std::cout<<"time after excuting while loop is "<<time4<<"\\n"<<std::endl;

}


//return outcome;

output.close();

}',depends="RcppArmadillo", plugins = "openmp")



