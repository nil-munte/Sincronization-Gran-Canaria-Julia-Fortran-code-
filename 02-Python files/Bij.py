with open('parameters_mesures.dat','r') as file:
    lines = file.readlines()

index = lines.index('L=39 ! Number of transmission lines\n')
L = lines[index].split()[0]
L = int(''.join(char for char in L if char.isdigit()))

index = lines.index('1 1 2 6050.00 26.297\n')

Zij = []
distances = []


for i in range(L):
    distances.append(float(lines[index + i].split()[4]))
    if (i in [0,5,33,34,35,36,37,38]):
        Zij.append(220.0**2/float(lines[index + i].split()[3]))
    else:
        Zij.append(66.0**2/float(lines[index + i].split()[3]))

res = []
Zij_new = []
distances_real = []

for i in range(L):
    distances_real.append(Zij[i] / 0.4)
    Zij_new.append(0.4 * distances[i] / L) if distances[i] != 0 else res.append(0)
    res.append(Zij[i] / distances[i]) if distances[i] != 0 else res.append(0)

for i in range(len(Zij_new)):
    print(distances[i])
    print(distances_real[i])


