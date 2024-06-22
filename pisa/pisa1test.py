with open('pisa2_summary.txt','r') as fp:
    all_result = fp.read().split('-----------------------------------------------------------')
    all_result.pop(0)
for result in all_result:
    result = result.strip().split('\n')
    if result[-1] == 'structure DB read fault':
        with open('pisa2_result_retry.txt','a') as fp:
            fp.write(result[0]+'\n')
    elif float(result[-1]) >= -9 and len(result[2]) == 6:
        with open('pisa2_result_highdletaG.txt','a') as fp:
            fp.write(result[0]+'\n') 
    elif float(result[-1]) < -9 and len(result[2]) == 6:
        with open('pisa2_result_lowdletaG.txt','a') as fp:
            fp.write(result[0]+'\n') 
    else:
        with open('pisa2_result_retry.txt','a') as fp:
            fp.write(result[0]+'\n')