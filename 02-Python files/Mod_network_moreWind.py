import copy

n = 0 # number of new nodes
N = 1 # Multiplication wind factor (xN)

# Change some parameters of the network
Bij = 1
node_to_connect = 2 # node to which the new nodes are connected

with open('parameters.dat','r') as file:
    lines = file.readlines()

# Modify N and L
lines[1] = f'N={29 + n} ! Number of nodes\n'
lines[3] = f'L={39 + n} ! Number of transmission lines\n'
lines[4] = f'Nassets={6 + n} ! Number of assets\n'

index = lines.index('dispFile=\'dispatch/dispatch_07-02-2018.dat\' ! dispatch file\n')
lines[index] = f'dispFile=\'dispatch/dispatch_07-02-2018_moreWind_n{n}_N{N}_Bij{Bij}.dat\' ! dispatch file\n'

index = lines.index('resFile=\'results/frequency_07-02-2018.dat\' ! results file\n')
lines[index] = f'resFile=\'results/frequency_07-02-2018_moreWind_n{n}_N{N}_Bij{Bij}.dat\' ! results file\n'
lines.insert(index + 1, f'resFile_power=\'results/power_07-02-2018_n{n}_N{N}_Bij{Bij}.dat\'\n')
lines.insert(index + 2, f'resFile_demand=\'results/demand_07-02-2018.dat\'\n')

# Find the index where the new line should be added
index = lines.index('39 2 29 20000.00\n') + 1

for i in range(n):
    new_line = f'{40 + i} {node_to_connect} {30 + i} {Bij}\n'
    lines.insert(index + i, new_line)

# Find the index where the new line should be added
index = lines.index('          6    23     wind\n') + 1

for i in range(n):
    new_line = f'          {7 + i}    {30 + i}     wind\n'
    lines.insert(index + i, new_line)


# Write the modified content back to the file
with open(f'parameters_moreWind_n{n}_N{N}_Bij{Bij}.dat', 'w') as file:
    file.writelines(lines)

with open('dispatch/dispatch_07-02-2018.dat','r') as file:
    lines = file.readlines()

lines_copy = copy.copy(lines)

i = 0

suma_loads = 0.0
suma_refs = 0.0
suma_assets = 0.0
power_to_provide = 0.0

for index, line in enumerate(lines):

    if line == 'loads:\n':
        suma_loads = 0.0
        for j in range(1,30):
            suma_loads+=float(lines[index+j].split()[1])

    if line == 'Pref:\n':
        for j in range(n):
            lines_copy.insert(index + i + j,f'{30 + j} 0.00\n')
        suma_refs = 0.0
        pref = []
        for j in range(1,7):
            suma_refs+=float(lines[index+j].split()[1])
            pref.append((index + j,float(lines[index+j].split()[1])))
        i+=n

    if line == '* end dispatch *\n':

        suma_assets = 0.0
        for j in range(1,4):
            suma_assets+=float(lines[index-j].split()[1])*N

        pref = sorted(pref, key=lambda x: x[1])
        subt = (N-1) * suma_assets/N
        for k, (a, b) in enumerate(pref):
            pref[k] = (a, b - subt) if b - subt > 0 else (a, 0)
            subt = max(0, subt - b)
            lines_copy[a + i] = f'{lines[a].split()[0]} {float(pref[k][1]):.2f} {lines[a].split()[2]}\n'

        for j in range(1,4):
            lines_copy[index + i - j] = f'{lines[index-j].split()[0]} {suma_assets/(3 + n)}\n'
    
        for j in range(n):
            lines_copy.insert(index + i + j,f'{7 + j} {suma_assets/(3 + n)}\n')

        i+=n
        
with open(f'dispatch/dispatch_07-02-2018_moreWind_n{n}_N{N}_Bij{Bij}.dat', 'w') as file:
    file.writelines(lines_copy)
