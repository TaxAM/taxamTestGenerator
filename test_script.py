from os import system
commands = [
    'rm -r *.tsv',
    'clear',
    'python taxamTestGenerator -n pool_esc_a -s A,B -nt 9,9,9,9,9,9,9 -pt 0 -nr 100 -nc 20 -pm 0.85 -tr 100 -tc 100 -cr 0.75 -cc 0.90 -mc 0.65 > output.txt',
    'rm -r pool_tests',
    'mv pool_esc_a pool_tests'
]

for key, command in enumerate(commands):
    try:
        system(command)
        print(f'{key} -> {command}')
    except Exception as e:
        print(e)