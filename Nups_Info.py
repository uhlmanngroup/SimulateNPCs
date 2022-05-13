#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 16:22:39 2021

@author: maria
"""
import numpy as np
import math

#Nup133
class SelectNup:
    
    def __init__(self, nup, term, model):
        self.nup = nup.upper()    
        self.term = term.upper()
        self.model = model.upper()                
        self.ref = self.rotUnitCoords(self.nup, self.term, self.model)        
        self.z = self.distz(self.ref)
        self.r = self.radii(self.ref)
        self.ringAngles = self.allang(self.ref)
        self.rotAng, self.rotOffset = self.rotang(self.z, self.ringAngles)
     
   
    def rotUnitCoords(self, nup, term, model):
        #MODEL 5A9Q
        #Nup43 N getcr #1.1/0,R,9,I:4@CA
        if (model == "5A9Q"):
            if (nup == "NUP43"):
                if (term == "N"):
                    authN_0 = np.array([1196.307, 578.280, 1007.303])
                    authN_R = np.array([1147.239, 748.554, 456.064])
                    authN_9 = np.array([1152.966, 676.288, 983.055])
                    authN_I = np.array([1188.192, 850.048, 420.6544])
                    ref = np.array([authN_0, authN_R, authN_9, authN_I])     
             # Nup43 C getcr #1.1/0,R,9,I:380@CA
                if (term == "C"):            
                    authC_0 = np.array([1186.279, 583.369, 1005.809])
                    authC_R = np.array([1137.054, 744.396, 458.835])
                    authC_9 = np.array([1142.872, 680.005, 979.448])
                    authC_I = np.array([1178.411, 844.515, 422.209])
                    ref = np.array([authC_0, authC_R, authC_9, authC_I])
        
            if (nup == "NUP160"):
                if (term == "N"):
            # Nup160 N getcr #1.1/1,S,a,J:41@CA
                    authN_1 = np.array([1288.803, 610.106, 903.284])
                    authN_S = np.array([1252.870, 725.826, 537.068])
                    authN_a = np.array([1260.478, 703.239, 891.066])
                    authN_J = np.array([1292.376, 823.203, 513.853])
                    ref = np.array([authN_1, authN_S, authN_a, authN_J])
                if (term == "C"):
                    # Nup160 C getcr #1.1/1,S,a,J:1201@CA
                    authC_1 = np.array([1191.789, 611.125, 983.108])
                    authC_S = np.array([1138.511, 708.726, 470.174])
                    authC_a = np.array([1146.753, 724.288, 960.323])
                    authC_J = np.array([1187.019, 805.358, 443.700])
                    ref = np.array([authC_1, authC_S, authC_a, authC_J])
    
    
            if (nup == "NUP37"):
                if(term == "N"):
                    # Nup37 N getcr #1.1/2,T,b,K:9@CA
                    authN_2 = np.array([1237.157, 598.816, 976.333])
                    authN_T = np.array([1179.777, 724.079, 484.290])
                    authN_b = np.array([1185.984, 706.003, 941.802])
                    authN_K = np.array([1226.971, 822.146, 451.789])
                    ref = np.array([authN_2, authN_T, authN_b, authN_K])
                    
                if(term == "C"):
                    # Nup37 C getcr #1.1/2,T,b,K:326@CA
                    authC_2 = np.array([1236.402, 588.357, 965.286])
                    authC_T = np.array([1181.011, 735.291, 494.523])
                    authC_b = np.array([1187.724, 694.322, 932.184])
                    authC_K = np.array([1225.765, 831.973, 463.361])
                    ref = np.array([authC_2, authC_T, authC_b, authC_K])
        
            if(nup == "SEC13"):
                if(term == "N"):
              #Sec13 N getcr #1.1/6,X,F,O:14@CA   
                    authN_6 = np.array([1140.803, 763.184, 987.984])
                    authN_X = np.array([1197.434, 664.774, 1018.729])
                    authN_F = np.array([1192.476, 762.622, 409.475])
                    authN_O = np.array([1138.455, 663.969, 443.912])
                    ref = np.array([authN_6, authN_X, authN_F, authN_O])               
                if(term == "C"):
                    authC_6 = np.array([1135.524, 781.493, 987.406])
                    authC_X = np.array([1193.033, 683.319, 1018.399])
                    authC_F = np.array([1188.354, 744.028, 410.265])
                    authC_O = np.array([1133.460, 645.602, 444.980])
                    ref = np.array([authC_6, authC_X, authC_F, authC_O])
        #Sec13 C getcr #1.1/6,X,F,O:304@CA
    
            if(nup == "SEH1"):
                if(term == "N"):
        #Seh1 N getcr #1.1/7,Y,G,P:1@CA
                    authN_7 = np.array([1135.180, 687.944, 1024.433])
                    authN_Y = np.array([1188.146, 589.668, 1051.665])
                    authN_G = np.array([1181.088, 836.930, 376.589])
                    authN_P = np.array([1131.479, 740.769, 412.990])
                    ref = np.array([authN_7, authN_Y, authN_G, authN_P])
                if(term == "C"):
        # Seh1 C getcr #1.1/7,Y,G,P:322@CA
                    authC_7 = np.array([1116.284, 703.552, 1035.049])
                    authC_Y = np.array([1173.701, 606.808, 1066.186])
                    authC_G = np.array([1167.666, 818.661, 362.467])
                    authC_P = np.array([1113.311, 726.060, 400.072])
                    ref = np.array([authC_7, authC_Y, authC_G, authC_P])
        
            if(nup == "NUP85"):
                if(term == "N"):   # Nup85 N getcr #1.1/Q,8,Z,H:9@CA
                    authN_Q = np.array([1104.998, 733.983, 412.939])
                    authN_8 = np.array([1108.652, 694.300, 1022.656])
                    authN_Z = np.array([1162.782, 599.234, 1055.191])
                    authN_H = np.array([1156.268, 826.013, 373.120])
                    ref = np.array([authN_Q, authN_8, authN_Z, authN_H])
    
                        
                if(term == "C"): # Nup85 C Atom getcrd #1.1/Q,8,Z,H:475@CA
                    authC_Q = np.array([1139.629, 717.241, 431.447])
                    authC_8 = np.array([1143.905, 709.705, 1004.162])
                    authC_Z = np.array([1195.665, 611.320, 1030.803])
                    authC_H = np.array([1189.437, 816.319, 398.179])
                    ref = np.array([authC_Q, authC_8, authC_Z, authC_H])
    
    # Nup155 N getcr #1.1/A,B:20@CA
            if(nup == "NUP155"):
                if(term == "N"):
                    authN_A = np.array([1126.810, 704.682, 588.647])
                    authN_B = np.array([1125.768, 710.075, 833.438])
                    ref = np.array([authN_A, authN_B])
    
    # Nup 155 C getcr #1.1/A,B:730@CA
                if(term == "C"):
                    authC_A = np.array([1109.824, 708.842, 587.688])
                    authC_B = np.array([1108.682, 706.273, 834.057])
                    ref = np.array([authC_A, authC_B])
               
            if nup == "NUP133": # Most N-terminal
                if (term == "N"):
                    auth3 = np.array([1145.276, 1042.219, 915.997]) #Atom #1.1/3:518@CA 
                    authU = np.array([1243.443, 946.836, 960.514]) #Atom #1.1/U:518@CA 
                    authC = np.array([1238.959, 480.349, 470.415]) #Atom #1.1/C:518@CA 
                    authL = np.array([1139.371, 391.992, 512.747]) #Atom #1.1/L:518@CA 
                    ref = np.array([auth3, authU, authC, authL]) # angles might be wrong
                
                if(term == "C"):    # Nup133 C getcr #1.1/3,U,C,L:1156@CA
                    authC_3 = np.array([1202.376, 949.301, 952.715])
                    authC_U = np.array([1270.591, 837.073, 984.803])
                    authC_C = np.array([1265.747, 584.752, 449.213])
                    authC_L = np.array([1201.188, 478.470, 477.475])
                    ref = np.array([authC_3, authC_U, authC_C, authC_L])
                    
            
            if nup == "NUP96":
                if(term == "N"): # N terminus getcr #1.1/5,W,E,N:277@CA 
                    authN_5 = np.array([1140.382, 748.182, 998.210])
                    authN_W = np.array([1196.973, 649.829, 1029.073])
                    authN_E = np.array([1191.410, 777.390, 398.923])
                    authN_N = np.array([1137.431, 678.795, 433.476])
                    ref = np.array([authN_5, authN_W, authN_E, authN_N])
    
                if(term == "C"):#Nup96 getcr getcr #1.1/5,W,E,N:751@CA 
                    auth5 = np.array([1177.224, 782.184, 975.888])    
                    authE = np.array([1229.633, 746.212, 422.979])
                    authW = np.array([1233.821, 681.917, 1004.085])
                    authN = np.array([1175.575, 645.661, 454.899])
                    ref = np.array([auth5, authE, authW, authN])
            
                
            #"select #1.1-end/M,D,V,4:150"
            if nup == "NUP107": 
                if(term == "N"):
                    authM = np.array([1134.432, 521.102, 440.365])
                    authD = np.array([1195.877, 619.745, 406.082])
                    authV = np.array([1199.033, 807.615, 1024.308])
                    auth4 = np.array([1134.991, 905.906, 994.000])
                    ref = np.array([authM, authD, authV, auth4]) 
                
                if(term == "C"): #getcr #1.1/M,D,V,4:924@CA
                    authC_4 = np.array([1162.799, 987.654, 965.130])
                    authC_V = np.array([1243.897, 883.264, 1003.383])
                    authC_D = np.array([1235.108, 543.326, 426.048])
                    authC_M = np.array([1162.346, 441.103, 460.547])
                    ref = np.array([authC_4, authC_V, authC_D, authC_M])
            
            
        # MODEL 7PEQ
        if model == "7PEQ":
            if nup == "NUP133":
                #getcr /?C:518@CA
                if(term == "N"):
                    authN_AC = np.array([541.999, 1346.763, 823.525])
                    authN_BC = np.array([638.117, 1440.660, 862.794])
                    authN_CC = np.array([523.924, 970.110, 1319.985])
                    authN_DC = np.array([591.148, 851.023, 1279.022])
                    ref = np.array([authN_AC, authN_BC, authN_CC, authN_DC])
                    
                if term == "C":
                    #getcrd /?C:1156@CA
                    authC_AC = np.array([519.158, 1237.688, 796.818])
                    authC_BC = np.array([594.420, 1338.527, 834.656])
                    authC_CC = np.array([521.047, 1080.295, 1351.747])
                    authC_DC = np.array([567.701, 958.772, 1310.625])
                    ref = np.array([authC_AC, authC_BC, authC_CC, authC_DC])
                    
            if nup == "NUP107":
                #getcrd /?D:150@CA
                if(term == "N"):
                    authN_AD = np.array([600.839, 1200.739, 781.033])
                    authN_BD = np.array([666.899, 1288.639, 811.276])
                    authN_CD = np.array([606.785, 1099.072, 1375.937])
                    authN_DD = np.array([645.965, 989.664, 1345.110])
                    ref = np.array([authN_AD, authN_BD, authN_CD, authN_DD])
                 
                #getcrd /?D:924@CAAtom 
                if(term == "C"):
                    authC_AD = np.array([551.473, 1282.074, 784.417])
                    authC_BD = np.array([633.087, 1377.111, 821.083])
                    authC_CD = np.array([542.705, 1029.732, 1363.691])
                    authC_DD = np.array([595.238, 911.805, 1324.919])
                    ref = np.array([authC_AD, authC_BD, authC_CD, authC_DD])
                    
            if nup == "NUP96":
                if(term == "N"):
                    #getcrd /?E:333@CA
                    authN_AE = np.array([581.129, 1089.071, 800.177])
                    authN_BE = np.array([631.393, 1180.644, 831.500])
                    authN_CE = np.array([620.007, 1211.549, 1351.841])
                    authN_DE = np.array([640.068, 1104.186, 1324.706])
                    ref = np.array([authN_AE, authN_BE, authN_CE, authN_DE])
                if term == "C":
                    #getcrd /?E:922@CA
                    authC_AE = np.array([556.616, 1016.471, 797.658])
                    authC_BE = np.array([592.129, 1113.714, 829.995])
                    authC_CE = np.array([596.265, 1286.142, 1356.496])
                    authC_DE = np.array([619.069, 1179.706, 1320.086])
                    ref = np.array([authC_AE, authC_BE, authC_CE, authC_DE])
                    
            if nup == "SEC13":
                if(term == "N"):
                    #getcrd /?F:14@CA
                    authN_AF = np.array([600.079, 1050.919, 795.351])
                    authN_BF = np.array([643.468, 1140.012, 825.169])
                    authN_CF = np.array([640.796, 1248.753, 1356.790])
                    authN_DF = np.array([662.173, 1140.940, 1325.852])
                    ref = np.array([authN_AF, authN_BF, authN_CF, authN_DF])
                if term == "C":
                    #getcrd /?F:304@CA
                    authC_AF = np.array([601.219, 1069.828, 795.637])
                    authC_BF = np.array([647.763, 1158.457, 825.743])
                    authC_CF = np.array([641.149, 1229.836, 1355.791])
                    authC_DF = np.array([661.812, 1122.034, 1327.058])
                    ref = np.array([authC_AF, authC_BF, authC_CF, authC_DF])
                    
            if nup == "SEH1":
                if term == "N":
                    #getcrd /?G:5@CA
                    authN_AG = np.array([627.178, 991.967, 749.919])
                    authN_BG = np.array([658.435, 1078.524, 777.520])
                    authN_CG = np.array([673.007, 1304.920, 1402.417])
                    authN_DG = np.array([695.554, 1200.844, 1365.491])
                    ref = np.array([authN_AG, authN_BG, authN_CG, authN_DG])
                
                if term == "C":
                    #getcrd /?G:320@CA
                    authC_AG = np.array([638.093, 1000.623, 749.846])
                    authC_BG = np.array([670.636, 1085.255, 777.118])
                    authC_CG = np.array([683.546, 1295.839, 1401.506])
                    authC_DG = np.array([705.751, 1191.343, 1365.910])
                    ref = np.array([authC_AG, authC_BG, authC_CG, authC_DG])     
                    
            if nup == "NUP85":
                if term == "N":
                    #getcrd /?H:20@CA
                    authN_AH = np.array([636.843, 1010.590, 740.489])
                    authN_BH = np.array([670.745, 1095.561, 768.025])
                    authN_CH = np.array([682.540, 1285.526, 1410.533])
                    authN_DH = np.array([704.151, 1182.196, 1376.040])
                    ref = np.array([authN_AH, authN_BH, authN_CH, authN_DH])
                    
                if term == "C":
                    #getcrd /?H:651@CA
                    authC_AH = np.array([570.694, 994.023, 825.221])
                    authC_BH = np.array([606.277, 1088.162, 855.288])
                    authC_CH = np.array([612.060, 1307.675, 1330.730])
                    authC_DH = np.array([636.318, 1197.414, 1292.414])
                    ref = np.array([authC_AH, authC_BH, authC_CH, authC_DH])
                    
            if nup == "NUP43":
                if term == "N":
                    #getcrd /?I:4@CA
                    authN_AI = np.array([607.400, 977.055, 801.079])
                    authN_BI = np.array([638.557, 1065.840, 829.237])
                    authN_CI = np.array([650.766, 1322.411, 1353.145])
                    authN_DI = np.array([675.070, 1213.263, 1313.972])
                    ref = np.array([authN_AI, authN_BI, authN_CI, authN_DI])
                    
                if term == "C":
                    #getcrd /?I:380@CA
                    authC_AI = np.array([615.244, 984.988, 802.829])
                    authC_BI = np.array([647.684, 1072.310, 830.776])
                    authC_CI = np.array([658.171, 1314.277, 1350.619])
                    authC_DI = np.array([682.206, 1204.628, 1312.594])
                    ref = np.array([authC_AI, authC_BI, authC_CI, authC_DI])
                    
            if nup == "NUP160":
                if term == "N":
                    #getcrd /?J:78@CA
                    authN_AJ = np.array([521.159, 984.775, 890.577])
                    authN_BJ = np.array([558.559, 1085.539, 922.612])
                    authN_CJ = np.array([558.979, 1321.315, 1268.981])
                    authN_DJ = np.array([545.739, 1234.093, 1248.102])
                    ref = np.array([authN_AJ, authN_BJ, authN_CJ, authN_DJ])
                    
                if term == "C":
                    #getcrd /?J:1195@CA
                    authC_AJ = np.array([590.711, 1029.030, 807.313])
                    authC_BJ = np.array([630.923, 1119.928, 837.190])
                    authC_CJ = np.array([632.207, 1271.078, 1346.311])
                    authC_DJ = np.array([654.672, 1162.184, 1312.629])
                    ref = np.array([authC_AJ, authC_BJ, authC_CJ, authC_DJ])
                    
            if nup == "NUP37":
                if term == "N":
                    #getcrd /?K:18@CA
                    authN_AK = np.array([533.337, 992.853, 833.045])
                    authN_BK = np.array([569.551, 1092.930, 864.748])
                    authN_CK = np.array([574.306, 1310.683, 1325.339])
                    authN_DK = np.array([590.786, 1213.351, 1280.732])
                    ref = np.array([authN_AK, authN_BK, authN_CK, authN_DK])
                    
                if term == "C":
                    #getcrd /?K:324@CA
                    authC_AK = np.array([553.058, 985.968, 839.464])
                    authC_BK = np.array([588.090, 1082.696, 870.158])
                    authC_CK = np.array([593.860, 1317.043, 1317.931])
                    authC_DK = np.array([600.402, 1212.340, 1261.130])
                    ref = np.array([authC_AK, authC_BK, authC_CK, authC_DK])
                    

        if model == "7PER":                    
            if nup == "NUP205":
                if term == "N":
                    #getcrd /D,J,V,P:9@CA
                    authN_D = np.array([766.579, 1044.574, 1052.258])
                    authN_J = np.array([708.820, 1125.978, 1045.733])
                    authN_V = np.array([686.322, 1092.394, 1186.438])
                    authN_P = np.array([753.571, 1153.951, 1195.119])
                    ref = np.array([authN_D, authN_J, authN_V, authN_P])
                                    
                if term == "C":       
                    #getcrd /D,J,V,P:1692@CA
                    authC_D = np.array([683.834, 1148.509, 1097.252])
                    authC_J = np.array([742.309, 1231.826, 1045.529])
                    authC_V = np.array([717.704, 986.583, 1198.460])
                    authC_P = np.array([670.166, 1058.951, 1134.355])
                    ref = np.array([authC_D, authC_J, authC_V, authC_P])
                    
            if nup == "NUP155":
                if term == "N":
                    #getcrd /E,K,Q,W:20@CA
                    authN_E = np.array([668.125, 986.789, 1118.243])
                    authN_K = np.array([649.617, 1074.233, 1045.931])
                    authN_Q = np.array([676.362, 1246.211, 1104.215])
                    authN_W = np.array([638.337, 1131.825, 1177.548])
                    ref = np.array([authN_E, authN_K, authN_Q, authN_W])
                if term == "C":                    
                    #getcrd /E,K,Q,W:1375@CA
                    authC_E = np.array([701.274, 1086.404, 1084.926])
                    authC_K = np.array([718.793, 1165.954, 1023.777])
                    authC_Q = np.array([708.040, 1150.287, 1154.007])
                    authC_W = np.array([698.384, 1045.048, 1205.549])
                    ref = np.array([authC_E, authC_K, authC_Q, authC_W])
                    
            if nup == "NUP93":
                if term == "N":
                    #getcrd /C,I,O,U:173@CA
                    authN_C = np.array([679.385, 1115.497, 1107.962])
                    authN_I = np.array([693.222, 1207.092, 1063.216])
                    authN_O = np.array([683.229, 1104.768, 1124.793])
                    authN_U = np.array([672.276, 1012.543, 1162.728])
                    ref = np.array([authN_C, authN_I, authN_O, authN_U])
                if term == "C":     
                    #getcrd /C,I,O,U:815@CA
                    authC_C = np.array([664.864, 1048.350, 1112.871])
                    authC_I = np.array([663.228, 1147.690, 1045.459])
                    authC_O = np.array([669.738, 1171.469, 1114.178])
                    authC_U = np.array([633.985, 1069.780, 1161.499])
                    ref = np.array([authC_C, authC_I, authC_O, authC_U])
                    
            if nup == "NUP54":
                if term == "N":
                    #getcrd /F,L,R,X:128@CA
                    authN_F = np.array([793.780, 1098.840, 1103.565])
                    authN_L = np.array([796.577, 1199.915, 1102.008])
                    authN_R = np.array([785.953, 1096.559, 1151.106])
                    authN_X = np.array([782.070, 996.679, 1138.960])
                    ref = np.array([authN_F, authN_L, authN_R, authN_X])
                    
                if term == "C":   
                    #getcrd /F,L,R,X:493@CA
                    authC_F = np.array([733.940, 1032.037, 1143.314])
                    authC_L = np.array([744.925, 1117.484, 1114.681])
                    authC_R = np.array([741.041, 1166.010, 1098.352])
                    authC_X = np.array([746.295, 1087.568, 1129.866])
                    ref = np.array([authC_F, authC_L, authC_R, authC_X])
                    
            if nup == "P58P45":
                if term == "N":
                    #getcrd /G,M,S,Y:248@CA
                    authN_G = np.array([796.774, 1096.951, 1088.269])
                    authN_M = np.array([805.350, 1200.433, 1088.999])
                    authN_S = np.array([787.615, 1099.730, 1166.393])
                    authN_Y = np.array([789.007, 993.510, 1152.683])
                    ref = np.array([authN_G, authN_M, authN_S, authN_Y])
                if term == "C":    
                    #getcrd /G,M,S,Y:418@CA
                    authC_G = np.array([741.218, 1036.293, 1147.267])
                    authC_M = np.array([748.939, 1122.796, 1121.192])
                    authC_S = np.array([747.976, 1160.386, 1095.709])
                    authC_Y = np.array([750.157, 1082.239, 1123.280])
                    ref = np.array([authC_G, authC_M, authC_S, authC_Y])
                    
            if nup == "NUP62":
                if term == "N":
                    #getcrd /H,N,T,Z:334@CA
                    authN_H = np.array([782.857, 1100.463, 1096.103])
                    authN_N = np.array([789.127, 1199.869, 1090.957])
                    authN_T = np.array([774.205, 1097.290, 1157.357])
                    authN_Z = np.array([773.444, 997.033, 1149.114])
                    ref = np.array([authN_H, authN_N, authN_T, authN_Z])
                if term == "C":       
                    #getcrd /H,N,T,Z:502@CA
                    authC_H = np.array([727.667, 1037.177, 1144.194])
                    authC_N = np.array([737.704, 1120.965, 1113.178])
                    authC_T = np.array([734.202, 1161.708, 1097.243])
                    authC_Z = np.array([738.479, 1085.264, 1130.228])
                    ref = np.array([authC_H, authC_N, authC_T, authC_Z])
                    
        if model == "5IJN":
            
            if nup == "NUP205":
                if term == "N":
                    #getcrd #1.1/D,J,P,V:9@CA
                    authN_D = np.array([960.284, 652.509, 782.339])
                    authN_J = np.array([961.280, 678.281, 811.489])
                    authN_P = np.array([960.962, 767.045, 644.009])
                    authN_V = np.array([964.721, 743.088, 614.105])
                    ref = np.array([authN_D, authN_J, authN_P, authN_V])    
                    
                if term == "C":
                    #getcrd #1.1/D,J,P,V:1692@CA
                    authC_D = np.array([1046.501, 746.245, 721.493])
                    authC_J = np.array([992.239, 816.461, 788.789])
                    authC_P = np.array([1038.535, 667.639, 704.642])
                    authC_V = np.array([990.301, 598.984, 634.525])
                    ref = np.array([authC_D, authC_J, authC_P, authC_V])
                    
            if nup == "NUP54":
                if term == "N":
                    #getcrd #1.1/F,L,R,X:128@CA
                    authN_F = np.array([934.656, 716.059, 736.701])
                    authN_L = np.array([938.484, 806.287, 742.363])
                    authN_R = np.array([932.467, 700.520, 695.076])
                    authN_X = np.array([935.376, 612.148, 686.866])
                    ref = np.array([authN_F, authN_L, authN_R, authN_X])   
                    
                if term == "C":
                    #getcrd #1.1/F,L,R,X:493@CA
                    authC_F = np.array([984.530, 653.126, 680.352])
                    authC_L = np.array([970.364, 714.857, 726.627])
                    authC_R = np.array([987.751, 762.904, 746.799])
                    authC_X = np.array([967.570, 703.607, 701.773])
                    ref = np.array([authC_F, authC_L, authC_R, authC_X])

                
            if nup == "NUP58":
                if term == "N":
                    #getcrd #1.1/G,M,S,Y:248@CA
                    authN_G = np.array([930.053, 710.449, 750.623])
                    authN_M = np.array([930.742, 808.497, 755.842])
                    authN_S = np.array([927.096, 706.032, 681.392])
                    authN_Y = np.array([928.412, 609.629, 673.023])
                    ref = np.array([authN_G, authN_M, authN_S, authN_Y])
                if term == "C":     
                    #getcrd #1.1/G,M,S,Y:418@CA
                    authC_G = np.array([978.520, 659.710, 677.663])
                    authC_M = np.array([967.247, 721.001, 720.364])
                    authC_S = np.array([981.755, 756.542, 750.005])
                    authC_Y = np.array([964.003, 697.693, 708.018])
                    ref = np.array([authC_G, authC_M, authC_S, authC_Y])
            
            if nup == "NUP155":
                if term == "N":
                    #getcrd #1.1/A,B,E,K,Q,W:20@CA
                    authN_A = np.array([1127.765, 714.280, 844.037])
                    authN_B = np.array([1129.475, 701.251, 582.893])
                    authN_E = np.array([1035.081, 602.042, 702.602])
                    authN_K = np.array([1091.232, 679.306, 777.270])
                    authN_Q = np.array([1037.108, 812.567, 723.198])
                    authN_W = np.array([1093.999, 718.479, 650.195])
                    ref = np.array([authN_A, authN_B, authN_E, authN_K, authN_Q, authN_W])
                    
                if term == "C": 
                    #getcrd #1.1/A,B:863@CA, getcrd #1.1/E,K,Q,W:1375@CA 
                    authC_A = np.array([1123.935, 733.451, 853.606])
                    authC_B = np.array([1126.827, 681.928, 573.229])
                    authC_E = np.array([1050.620, 694.556, 764.925])
                    authC_K = np.array([1033.018, 765.870, 830.250])
                    authC_Q = np.array([1051.867, 722.745, 661.181])
                    authC_W = np.array([1031.677, 647.801, 595.052])
                    ref = np.array([authC_A, authC_B, authC_E, authC_K, authC_Q, authC_W])
            
            if nup == "NUP93":
                if term == "N":
                    #getcrd #1.1/C,I,O,U:1@CA
                    authN_C = np.array([976.397, 666.536, 729.429])
                    authN_I = np.array([961.761, 749.355, 763.917])
                    authN_O = np.array([975.962, 748.926, 698.630])
                    authN_U = np.array([960.689, 668.109, 665.067])
                    ref = np.array([authN_C, authN_I, authN_O, authN_U])
                    
                if term == "C": 
                    #getcrd #1.1/C,I,O,U:815@CA
                    authC_C = np.array([1060.023, 674.793, 712.464])
                    authC_I = np.array([1055.745, 758.808, 794.813])
                    authC_O = np.array([1057.252, 740.183, 714.214])
                    authC_U = np.array([1054.300, 653.797, 630.882])
                    ref = np.array([authC_C, authC_I, authC_O, authC_U])

            if nup == "NUP62":
                if term == "N":
                    #getcrd #1.1/H,N,T,Z:334@CA
                    authN_H = np.array([945.004, 713.228, 744.610])
                    authN_N = np.array([946.332, 804.477, 752.982])
                    authN_T = np.array([942.333, 702.920, 686.442])
                    authN_Z = np.array([943.862, 613.538, 676.683])
                    ref = np.array([authN_H, authN_N, authN_T, authN_Z])
                    
                if term == "C": 
                    #getcrd #1.1/H,N,T,Z:502@CA
                    authC_H = np.array([991.687, 656.938, 681.239])
                    authC_N = np.array([978.236, 716.679, 727.738])
                    authC_T = np.array([994.722, 758.872, 745.506])
                    authC_Z = np.array([975.468, 701.654, 701.185])
                    ref = np.array([authC_H, authC_N, authC_T, authC_Z])

                                     
        if model == "5IJO":                    
            if nup == "NUP54":
                if term == "N":
                    #getcrd /F,L,R,X:128@CA
                    authN_F = np.array([934.656, 716.059, 736.701])
                    authN_L = np.array([938.484, 806.287, 742.363])
                    authN_R = np.array([932.467, 700.520, 695.076])
                    authN_X = np.array([935.376, 612.148, 686.866])
                    ref = np.array([authN_F, authN_L, authN_R, authN_X])

                if term == "C": 
                    #getcrd /F,L,R,X:493@CA
                    authC_F = np.array([984.530, 653.126, 680.352])
                    authC_L = np.array([970.364, 714.857, 726.627])
                    authC_R = np.array([987.751, 762.904, 746.799])
                    authC_X = np.array([967.570, 703.607, 701.773])
                    ref = np.array([authC_F, authC_L, authC_R, authC_X])

            if nup == "NUP188":
                if term == "N":
                    #getcrd /J,V:1@CA
                    authN_J = np.array([963.691, 696.108, 808.948])
                    authN_V = np.array([964.642, 722.851, 612.857])
                    ref = np.array([authN_J, authN_V])
                    
                if term == "C": 
                    #getcrd /J,V:1564@CA
                    authC_J = np.array([989.042, 805.579, 797.104])
                    authC_V = np.array([987.194, 611.748, 624.551])
                    ref = np.array([authC_J, authC_V])
                    
            if nup == "NUP205":
                if term == "N":
                    #getcrd /D,P:9@CA
                    authN_D = np.array([960.284, 652.509, 782.339])
                    authN_P = np.array([960.962, 767.045, 644.009])
                    ref = np.array([authN_D, authN_P])
                    
                if term == "C": 
                    #getcrd /D,P:1692@CA
                    authC_D = np.array([1046.501, 746.245, 721.493])
                    authC_P = np.array([1038.535, 667.639, 704.642])
                    ref = np.array([authC_D, authC_P])

            if nup == "NUP155":
                if term == "N":
                    #getcrd /A,B,E,K,Q,W:20@CA
                    authN_A = np.array([1127.765, 714.280, 844.037]) # as in 5IJN 1.1
                    authN_B = np.array([1129.475, 701.251, 582.893])
                    authN_E = np.array([1035.081, 602.042, 702.602])
                    authN_K = np.array([1091.232, 679.306, 777.270])
                    authN_Q = np.array([1037.108, 812.567, 723.198])
                    authN_W = np.array([1093.999, 718.479, 650.195])
                    ref = np.array([authN_A, authN_B, authN_E, authN_K, authN_Q, authN_W])
                    
                if term == "C": 
                    #getcrd /A,B:863@CA, getcrd /E,K,Q,W:1375@CA
                    authC_A = np.array([1123.935, 733.451, 853.606])
                    authC_B = np.array([1126.827, 681.928, 573.229])
                    authC_E = np.array([1050.620, 694.556, 764.925])
                    authC_K = np.array([1033.018, 765.870, 830.250])
                    authC_Q = np.array([1051.867, 722.745, 661.181])
                    authC_W = np.array([1031.677, 647.801, 595.052])
                    ref = np.array([authC_A, authC_B, authC_E, authC_K, authC_Q, authC_W])
                    
            if nup == "NUP93":
                if term == "N":
                    #getcrd /C,I,O,U:1@CA
                    authN_C = np.array([976.397, 666.536, 729.429])
                    authN_I = np.array([961.761, 749.355, 763.917])
                    authN_O = np.array([975.962, 748.926, 698.630])
                    authN_U = np.array([960.689, 668.109, 665.067])
                    ref = np.array([authN_C, authN_I, authN_O, authN_U])
                    
                if term == "C": 
                    #getcrd /C,I,O,U:815@CA
                    authC_C = np.array([1060.023, 674.793, 712.464])
                    authC_I = np.array([1055.745, 758.808, 794.813])
                    authC_O = np.array([1057.252, 740.183, 714.214])
                    authC_U = np.array([1054.300, 653.797, 630.882])
                    ref = np.array([authC_C, authC_I, authC_O, authC_U])
                    
            if nup == "P58P45":
                if term == "N":
                    #getcrd /G,M,S,Y:248@CA
                    authN_G = np.array([930.053, 710.449, 750.623])
                    authN_M = np.array([930.742, 808.497, 755.842])
                    authN_S = np.array([927.096, 706.032, 681.392])
                    authN_Y = np.array([928.412, 609.629, 673.023])
                    ref = np.array([authN_G, authN_M, authN_S, authN_Y])
                    
                if term == "C":     
                    #getcrd /G,M,S,Y:418@CA
                    authC_G = np.array([978.520, 659.710, 677.663])
                    authC_M = np.array([967.247, 721.001, 720.364])
                    authC_S = np.array([981.755, 756.542, 750.005])
                    authC_Y = np.array([964.003, 697.693, 708.018])
                    ref = np.array([authC_G, authC_M, authC_S, authC_Y])
                
            if nup == "NUP62":
                if term == "N":
                    #getcrd /H,N,T,Z:334@CA
                    authN_H = np.array([945.004, 713.228, 744.610])
                    authN_N = np.array([946.332, 804.477, 752.982])
                    authN_T = np.array([942.333, 702.920, 686.442])
                    authN_Z = np.array([943.862, 613.538, 676.683])
                    ref = np.array([authN_H, authN_N, authN_T, authN_Z])
                    
                if term == "C": 
                    #getcrd /H,N,T,Z:502@CA
                    authC_H = np.array([991.687, 656.938, 681.239])
                    authC_N = np.array([978.236, 716.679, 727.738])
                    authC_T = np.array([994.722, 758.872, 745.506])
                    authC_Z = np.array([975.468, 701.654, 701.185])
                    ref = np.array([authC_H, authC_N, authC_T, authC_Z])
                    
        return ref 
    
    def centre(self, model):
        if model == "5A9Q":  
            c = np.array([711.36, 711.36]) # central axis 
            
        if model == "7PEQ" or model == "7PER": 
            c = np.array([1094.4, 1094.4]) # central axis 7PEQ determined from EM density map
            
        if model == "5IJN" or model == "5IJO" :
            c = np.array([716.9, 713.55])
            
        return c

    
    def cornang(self, p_ref, p1): #TODO: Does this also work for rotational units in a different quadrant? 
        "Angle between node p0, p1 (PDB 5A9Q, coordinates from chimera) and the centre"
        c = self.centre(self.model)
            
        p_refnew = p_ref[:2] - c   # change coordinate system to 0 centre
        p1new = p1[:2] - c
        
        rotateby = np.arctan2(p_refnew[1], p_refnew[0]) # rotate such that ref overlays with positive x axis 
        arctan2_p1 = np.arctan2(p1new[1], p1new[0])
        ang = arctan2_p1 - rotateby

        if ang < -math.pi:
            ang += 2*math.pi
        if ang > math.pi:
            ang -= 2*math.pi

#  #      if (abs(ang) > math.pi): ang = 2*math.pi - ang

        return ang#np.arctan2(p1new[1], p1new[0]) - rotateby

    
    def allang(self, ref):

        ringAngles = []
        for i in range(len(ref)):
            ringAngles.append(self.cornang(ref[0], ref[i]))
            
        ringAngles = ringAngles - np.mean(ringAngles)
        return ringAngles # note: If the rotation is counterclockwise, the angle has a positive measure
    
    def distz(self, p):
        
        z = np.zeros(len(p))
    
        for i in range(len(z)):
            z[i] = p[i][2]
    
        z -= min(z) # offset so that smallest z is 0 TODO: this needs to be global for multicolour 
        z /= 10 # Angstrom to nm 
        z = [round(i, ndigits = 2) for i in z] # round to nearest .1 nm
        return np.array(z)
    
    def radii(self, p):
        c = self.centre(self.model)

        r = np.zeros(len(p))
        for i in range(len(r)):
            x = p[i][0] - c[0]
            y = p[i][1] - c[1]
            r[i] = np.sqrt(x**2 + y**2)
    
            
        r /= 10
        r = [round(i, ndigits=1) for i in r]
        return r
    
    def rotang(self, z, ringAng):
        "Find azimuthal angle between nearest 'corners' of two rings"
        midplane = np.mean(z) 
        octOffset = 2*np.pi/8 # TODO: different sym?
        
        #Angle of rotational unit in CR. 
        crAng = np.mean(ringAng[z > midplane]) 
        crAngCW = crAng - octOffset # next rot unit CW
        crAngACW = crAng + octOffset # next rot unit ACW
        
        #Angle of rotational unit in NR
        nrAng = np.mean(ringAng[z < midplane])

        ang = [(nrAng - i) for i in [crAngCW, crAng, crAngACW]]

        minval, minindex = min([(abs(val), idx) for (idx, val) in enumerate(ang)])
        
        offset = [-octOffset, 0, octOffset][minindex]
        return ang[minindex], offset
        
        
