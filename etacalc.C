void etacalc(){
double theta, eta, dz,z_center,z_edge, theta_edge;
bool dz_computation;

static const double c_light     = 0.299792458 ; // mm/ps , or m/ns
static const double inv_c_light = 3.335640952 ; // ps/mm or ns/m

cout<<"insert eta:"<<endl;
cin>>eta;
theta = 2*TMath::ATan(exp(- eta)); // in rad
//theta = theta * 360 / 6.28;
cout<<"theta = "<<theta * 360 / 6.28<<endl;
cout<<"do you want to compute dz?  (1 for yes, 0 for no)"<<endl;
cin>>dz_computation;
if(dz_computation)  {
    //z_center = 1.29*sin(-theta);
    z_center = 1.29 * TMath::Tan ( TMath::Pi()/2. - theta);
    //z_edge = 1.29*sin(-theta - 2*TMath::ATan(exp(- 0.017)));
    theta_edge = 2*TMath::ATan(exp(- eta - 0.0085)),
    z_edge = 1.290 * TMath::Tan ( TMath::Pi()/2. - theta_edge);// + TMath::Pi()/2. - 2*TMath::ATan(exp(- 0.017 - eta)));
    dz = (z_edge - z_center);
    cout<<"dz = "<<dz<<"m"<<endl;
    cout<<"dz_center : "<<z_center<<endl;
    cout<<"dz_edge : "<< z_edge<<endl;
    cout<<"dt = "<<dz*1000*inv_c_light<<"ps"<<endl;
    }
}
