import sys
with open('results.txt', 'r') as f:
    for line in f.readlines():
        if line.startswith('Final:'):
            sizes = line.split()[3:7]
            break
    print(sizes)
    for s in sizes:
        if abs(float(s) - 14.65) > 0.2:
            print(f'{s} out of range!')
            sys.exit(1)

print('Uniform Sphere check - OK')
