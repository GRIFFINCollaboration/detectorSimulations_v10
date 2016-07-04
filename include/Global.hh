#ifndef GLOBAL_HH
#define GLOBAL_HH

inline G4double TransX(G4double x, G4double y, G4double z, G4double theta, G4double phi) {
    return (pow(x*x+y*y+z*z,0.5)*sin(theta)*cos(phi));
}

inline G4double TransY(G4double x, G4double y, G4double z, G4double theta, G4double phi) {
    return (pow(x*x+y*y+z*z,0.5)*sin(theta)*sin(phi));
}

inline G4double TransZ(G4double x, G4double y, G4double z, G4double theta) {
    return (pow(x*x+y*y+z*z,0.5)*cos(theta));
}

#endif
