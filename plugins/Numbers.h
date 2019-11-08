#ifndef NUMBERS_H
#define NUMBERS_H

// --------- For Kinematic Fit
inline float GetErrEt(float Et, float Eta) {
  float InvPerr2, a, b, c;
  if(fabs(Eta) < 1.4) {
    a = 0.00235;
    b = 0.796;
    c = 0.00161;
  }
  else {
    a = 0.0068;
    b = 1.07;
    c = 0.000905;
  }
  InvPerr2 = Et*Et*(a/(Et * Et) + b/Et + c);
  return InvPerr2;
}

inline float GetErrEta(float Et, float Eta) {
  float InvPerr2, a, b, c;
  if(fabs(Eta) < 1.4) {
    a = 0.789;
    b = 0.0167;
    c = 0;
  }
  else {
    a = 1.41;
    b = 0.0199;
    c = 0.;
  }
  InvPerr2 = a/(Et * Et) + b/Et + c;
  return InvPerr2;
}

inline float GetErrPhi(float Et, float Eta) {
  float InvPerr2, a, b, c;
  if(fabs(Eta) < 1.4) {
    a = 0.709;
    b = 0.0285;
    c = 0;
  }
  else {
    a = 0.679;
    b = 0.0269;
    c = 0;
  }
  InvPerr2 = a/(Et * Et) + b/Et + c;
  return InvPerr2;
}



#endif
