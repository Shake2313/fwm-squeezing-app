"""Gold-frozen constant coupling compared with variable Gaussian closure."""
import json
from pathlib import Path
import time
from .heldout import evaluate


def main():
    started=time.perf_counter()
    directory=Path(__file__).parent
    gaussian=json.loads((directory/'heldout_results.json').read_text(encoding='utf-8'))
    factor=.5594938027
    data={'constant_factor':factor,'fit_target':15.5,
          'role':'Single Gold-calibrated offdiagonal effective mixing participation; no separately identified Zeeman/spatial/Raman/density cause.',
          'cases':{},'nearby_sweeps':{}}
    for name,row in gaussian['cases'].items():
        item=evaluate(row['input'],constant_factor=factor)
        item['G_s_gaussian']=row['G_s_candidate']
        item['G_c_gaussian']=row['G_c_candidate']
        item['literature_targets']=row['literature_targets']
        data['cases'][name]=item
    for parameter,rows in gaussian['nearby_sweeps'].items():
        data['nearby_sweeps'][parameter]=[]
        for row in rows:
            item=evaluate(row['input'],constant_factor=factor)
            item['G_s_gaussian']=row['G_s_candidate']
            item['G_c_gaussian']=row['G_c_candidate']
            data['nearby_sweeps'][parameter].append(item)
    data['elapsed_seconds']=time.perf_counter()-started
    (directory/'alternative_results.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({
        'cases':{name:{k:v for k,v in row.items() if k!='input'} for name,row in data['cases'].items()},
        'sweeps':{name:[{'value':row['input'][name],'legacy':row['G_s_before'],
                        'constant':row['G_s_candidate'],'gaussian':row['G_s_gaussian']}
                        for row in rows] for name,rows in data['nearby_sweeps'].items()},
        'elapsed_seconds':data['elapsed_seconds']},indent=2))


if __name__ == '__main__':
    main()
