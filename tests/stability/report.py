#! /usr/bin/python3

'''
CQ_Gears - CadQuery based involute profile gear generator

Copyright 2021 meadiode@github

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
'''

import os
import json
from datetime import datetime
from collections import Counter, defaultdict

FAIL_TAGS = ('UNDEFINED', 'PRECALC', 'BUILDING', 'NOT_SOLID', 'VOLUME',
             'BBOX_Z', 'BBOX_XY', 'TIMEOUT')


class Report:
    '''Summarizes test result data.

    Also organizes parameters which failed to produce a valid gear for further
    case reproduction and debugging.
    '''
    
    def __init__(self, filename):
        info = filename.split('_')
        self.branch = info[-3]
        self.git_sha = info[-2]
        date_str = info[-1].split('.')[0]        
        self.dt = datetime.strptime(date_str, '%Y%m%d%H%M%S')
        
        data = json.loads(open(filename).read())
        
        self.n_passed = data['summary']['passed']
        self.n_failed = data['summary']['failed']
        self.n_total = data['summary']['total']
        self.duration = data['duration']
        
        self.fail_tags = Counter()
        self.fail_tags_per_test = defaultdict(Counter)
        self.fails = defaultdict(lambda: defaultdict(list))
        self.perc_per_test = defaultdict(Counter)
        
        for test in data['tests']:
            test_name = test['nodeid'].split('::')[-2]

            self.perc_per_test[test_name]['total'] += 1

            if test['outcome'] == 'passed':
                continue

            self.perc_per_test[test_name]['failed'] += 1

            tag = test['metadata']['failure_tag']
            self.fail_tags[tag] += 1
            
            params = test['metadata']['gear_params']
            
            self.fails[test_name][tag].append(params)
            self.fail_tags_per_test[test_name][tag] += 1
            
            
    def print_summary(self):
        dt = self.dt.strftime('%d-%m-%Y, %H:%M:%S')
        prc = self.n_passed / self.n_total * 100.0
        fprc = 100.0 - prc
        
        secs = int(self.duration)
        hrs = secs // (60 * 60)
        secs %= (60 * 60)
        mins = secs // 60
        secs %= 60
        
        print(f'branch:        {self.branch}')
        print(f'git SHA:       {self.git_sha}')
        print(f'date & time:   {dt}')
        print(f'running time:  {hrs} hours, {mins} minutes, {secs} seconds')
        print(f'passed:        {self.n_passed} ({prc:.2f}%)')
        print(f'failed:        {self.n_failed} ({fprc:.2f}%)')
        
        header = ' ' * 26
        for tag in FAIL_TAGS:
            header = header + f'{tag:>10}'

        header = header + '  % FAILED'
        
        print(header)
                
        for tname, tags in self.fail_tags_per_test.items():
            line = f'{tname:<26}'
            for tag in FAIL_TAGS:
                n = tags[tag]
                line = line + f'{n:>10}'

            tot = self.perc_per_test[tname]['total']
            tot_fails = self.perc_per_test[tname]['failed']
            perc = tot_fails / tot * 100.0
    
            line = line + f'{perc:>10.1f}'

            print(line)


    @staticmethod
    def gather_reports(dir_):
        reps = []

        for fname in os.listdir(dir_):
            fname = os.path.join(dir_, fname)
            if os.path.isfile(fname) and fname.endswith('.json'):
                reps.append(Report(fname))
                
        reps.sort(key=lambda e: e.dt, reverse=True)

        return reps


if __name__ == '__main__':
    reps = Report.gather_reports('./reports')

    print('\n')
    for rep in reps[::-1]:
        rep.print_summary()
        print()
        print('=' * 116)
        print('=' * 116)
        print()
