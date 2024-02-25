# scaling according to Bosia2013
import numpy as np
import math
import matplotlib.pyplot as plt

#Plot setting:
plt.rcParams.update({'font.size': 8})
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman"]

#load SRIM Output
File = open("Srim\\VACANCY_Z2_1AN.txt", "r")
Input = File.read().splitlines()
File.close()

#delete file header
a = 0
while "---------- " not in Input[a]:
    a = a + 1
del (Input[0:a + 1])
del (Input[(len(Input) - 2):len(Input)])

#convert data to float
Data = []
for i in range(len(Input)):
    line = []
    j = len(Input[i].split())
    line.append(float(Input[i].split()[0].replace(",", ".")))
    line.append(float(Input[i].split()[j - 1].replace(",", ".")))
    Data.append(line)

Data2 = np.array(Data)
Faktor = np.array([[0.1, 0], [0, 1e8]])
Data3 = np.dot(Data2, Faktor)  # conversion in x [nm] and Vac/(cm*Ion)

F = 2.25e17  # Dose pro cm^2, to be edited as needed
D_vac_crit = 2.8e22  # critical vacancy density [vac/cm^3], Fairchild2012 uncertainty: +-0,095e22
Alpha = 4.4e22  # defect recombination probability factor [1/cm^3] from Bosia2013
Roh_D = 3.52  # diamond density [g/cm^3]
Roh_A=2.14 #minimum density [g/cm^3] measurement Bosia2013 - 2.1 according to Fairchild2012
Dz = Data3[2, 0] - Data3[1, 0]  # step width [nm]

n_ion=F/(Data3[-1,0]*10**(-7)) #ion/cm^3
n_atom=Roh_D/(12.11*1.661*10**(-24)) #Atom number
ppm_ion=n_ion*10**6/n_atom #in ppm
print('ppm_ion',"{:e}".format(ppm_ion))

#density calculation
Data4 = np.array([[0, 0, 0, 0, 3.52, 0]])
for i in range(len(Data3)):
    line = Data3[i, :]
    D_orig = Data3[i, 1] * F  # Vaccancy density from SRIM
    D_sat = Alpha * (1 - math.exp(-Data3[i, 1] * F / Alpha))  # Vaccancy density from Bosia2013
    Roh_F = Roh_D - (Roh_D - Roh_A) * (1 - math.exp(-Data3[i, 1] * F / Alpha))
    line2 = np.append(line, D_orig)
    line3 = np.append(line2, D_sat)
    line4 = np.append(line3, Roh_F)
    Z_new = Data4[i, 5] + Roh_D / Roh_F * Dz
    line5 = np.append(line4, Z_new)
    Data4 = np.append(Data4, [line5], axis=0)

Damage=0
i=0
pre=0
j=0
istop=Data4[-1,0]
istart=Data4[0, 2]
step=(Data4[-1, 5]-Data4[0, 5])/len(Data4)

#swelling calculation:
for i in range(len(Data4)):
    if pre<(D_vac_crit) and Data4[i,3] > (D_vac_crit):
        istart=Data4[i, 5]
        print(istart,'start')
    if pre>(D_vac_crit) and Data4[i,3] < (D_vac_crit):
    #if Data4[i,2] == 0 and pre > 0:
        istop=Data4[i, 5]
        if i == 0:
            istop = Data4[-1,5]
        print(istop,'stop')
    if Data4[i,3] > (D_vac_crit):
        Damage=Damage+Data4[i,3]
        j=j+1

    pre=Data4[i,3]

mean_Damage=Damage/j*(istop-istart)*10**(-7) #vaccancy concentration

#damagesquare=Damage*((istop)*10**(-7))
print('meanDamage in graphit layer v/cm²=',"{:e}".format(mean_Damage))
print('damage soll v/cm²',F*4)
#print('graphit thickness=',(istop-istart))
swelling = Data4[-1, 5] - Data4[-1, 0]
print('swelling=',swelling)

#save data
np.savetxt("VACANCY_rescaled2.txt", Data4, delimiter=";", header="z_orgi,Vac_orgi,Vac_F_orgi,Vac_F_sat,Roh_F,z_new")


# plot figure
fig, ax = plt.subplots(figsize=(4.0, 3))
plt.plot(Data4[:, 5], Data4[:, 3], label="Dose scaled model")
plt.plot(Data4[:, 0], Data4[:, 2], "k:", label="SRIM data")
plt.plot([0, max(Data4[:, 5])], [D_vac_crit+0.095e22, D_vac_crit+0.095e22],'gray',alpha=0.01)
plt.plot([0, max(Data4[:, 5])], [D_vac_crit-0.095e22, D_vac_crit-0.095e22],'gray',alpha=0.01)
ax.axhspan(D_vac_crit-0.095e22, D_vac_crit+0.095e22, xmin=0, xmax=500, alpha=0.5, color='C0',label="Graphitization threshold")
plt.ylabel("Vacancy concentration [1/cm³]")
plt.xlabel("Depth [nm]")
plt.legend()
plt.savefig('vacsim.pdf', transparent=True,bbox_inches ="tight")
plt.show()

