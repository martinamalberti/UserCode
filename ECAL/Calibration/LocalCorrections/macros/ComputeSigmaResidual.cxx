double ComputeSigmaResidual(double dm1, double sigma1, double dm2, double sigma2)
{

  double res = 0.5*(91.188+dm1+91.188+dm2) * sqrt( pow(sigma1/(dm1+91.188),2) - pow(sigma2/(dm2+91.188),2)  );
  cout << pow(sigma1/(dm1+91.188),2) << "  " << pow(sigma2/(dm2+91.188),2) << endl;
  return (res);
}
