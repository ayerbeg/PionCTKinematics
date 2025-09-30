# PionCTKinematics

This is a ROOT script made with ChatGPT-5 (-4 was not even able to understand the request). After several iteractions and corrections made by hand, the script is able to reproduce, quite exactly, the kinematics table of the Pion CT experiment in Hall C. The purpose is not only reproduce, but to find other values which allow to modify the table according different values of the Beam energy, Q2 and/or t. 

This needs to be validated yet. 
```
 Usage examples (from a shell):
    root -l -q "pionCT_kinematics.C(Q2,Ebeam,t)"
for example:
    root -l -q "pionCT_kinematics.C(5.0,11.0,-0.40)"
    root -l
    .L pionCT_kinematics.C+
    pionCT_kinematics(5.0,11.0,-0.40);
