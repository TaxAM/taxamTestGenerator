from os import system

commands = [
    'clear',
    'python taxamTestGenerator -n pool_esc_a -s A,B -nt 9,9,9,9,9,9,9 -pt 0 -nr 1000 -nc 1000 -pm 0.85 -tr 3000 -tc 1000 -cr 0.75 -cc 0.90 -mc 0.65 > output.txt',
    'rm -r pool_tests',
    'mv pool_esc_a pool_tests'
]

for command in commands:
    try:
        system(command)
    except Exception as e:
        print(e)