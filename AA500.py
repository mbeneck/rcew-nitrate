import pandas as pd
import matplotlib.pyplot as plt

class AA500_Result:
    def __init__(self, result_path, samplelist_path, master_path, result_mapping={'Results 1':'Nitrate', 'Results 2':'Phosphate', 'Results 3':'Ammonium'}, isv_thresholds = {'Nitrate':.005, 'Phosphate':.004, 'Ammonium': .005}, bbv_thresholds = {'Nitrate':.1, 'Phosphate':.1, 'Ammonium': .1}, pH_threshold = 8, spike = .1):
        self._isv_thresholds = isv_thresholds
        self._bbv_thresholds = bbv_thresholds
        self._pH_threshold = pH_threshold

        self._result_mapping = result_mapping
        self._raw_result_df = self._read_AA500_results(result_path)
        self._metadata = self._read_AA500_metadata(result_path)
        self._samplelist_df = self._read_samplelist_df(samplelist_path, master_path)
        self.result_df = self._merge_results_samplelist(result_mapping)
        self.result_df['AA500 Run Date'] = pd.to_datetime(self._metadata.loc['DATE', 'Value']+ ' ' + self._metadata.loc['TIME', 'Value'])
        self.result_df['AA500 Operator'] = self._metadata.loc['OPER', 'Value']
        self._calc_in_sample_std_QA()
        self._check_pH()
        self._calc_bbv()
        self.unspiked_result_df = self.result_df.loc[(slice(None), 0),:]
        if 1 in self.result_df.index.get_level_values('Spike'):
            self.spiked_result_df = self.result_df.loc[(slice(None), 1), :]
        else:
            self.spiked_result_df = pd.DataFrame()
        self._calc_recovery(spike)

    @staticmethod
    def _color_recovery(val):
        if 90 <= val <= 110:
            return 'background-color: green'
        elif 80 <= val <= 120:
            return 'background-color: orange'
        else:
            return 'background-color: red'

    def _calc_recovery(self, spike):
        if self.spiked_result_df.empty:
            self.recovery = pd.DataFrame()
        else:
            columns = [val + ' mean' for val in list(self._result_mapping.values())]
            unspiked = self.unspiked_result_df.reset_index().set_index('Sample ID')[columns]
            spiked = self.spiked_result_df.reset_index().set_index('Sample ID')[columns]
            recovery = (spiked - unspiked)/spike*100
            recovery = recovery.dropna()
            self.recovery = recovery
            self.recovery_styled = recovery.style.map(self._color_recovery)


    def _read_AA500_metadata(self, path):
        metadata = pd.read_csv(path, nrows=7, names=['Param', 'Value'], usecols=[0,1]).set_index('Param')
        return metadata

    def _read_AA500_results(self, path):
        df = pd.read_csv(path, header=14, parse_dates=['Date Time Stamp']).rename(columns = self._result_mapping)
        return df

    def _read_samplelist_df(self, samplelist_path, master_path):
        sample_df = pd.read_excel(samplelist_path)
        master_df = pd.read_excel(master_path)
        master_df['Sample Datetime'] = pd.to_datetime(master_df['Sample Collected'].astype(str) + ' ' + master_df['Time Sample Collected'].fillna('00:00:00').astype(str), errors='coerce')
        df=  sample_df.set_index('Sample ID').join(master_df.set_index('Sample ID').drop(columns='Cup Number'), how='left').reset_index().set_index('Cup Number')
        return df
    
    def _merge_results_samplelist(self, result_mapping):
        samples = self._raw_result_df[(self._raw_result_df['Cup Type'] == 'DSAMP')| (self._raw_result_df['Cup Type'] == '3SAMP') | (self._raw_result_df['Cup Type'] == 'SAMP')] 
        sample_results = samples.groupby('Cup Number')[list(result_mapping.values())].describe()
        sample_results= sample_results.loc[:, pd.IndexSlice[:, ['mean', 'std']]]

        err = sample_results.loc[:, pd.IndexSlice[:, 'std']].apply(lambda x: 2 * x)         ## add err (2* std)
        err.columns = pd.MultiIndex.from_tuples([(col[0], 'err') for col in err.columns])
        sample_results = pd.concat([sample_results, err], axis=1).sort_index(axis=1)

        sample_results.columns = sample_results.columns.map(' '.join)
        return self._samplelist_df.join(sample_results, on = 'Cup Number', how='right', validate='one_to_one').reset_index().set_index(['Sample ID', 'Spike']).sort_values('Sample ID')
    
    def _concat_qa_strings(self, df, value, column2):
        result = pd.Series()
        if value + ' QA' in df.columns:
            column1  = df[value + ' QA']
            result = pd.Series([
                ','.join([a for a in [x, y] if a])
                for x, y in zip(column1, column2)
            ], index = df.index)
        else:
            result = column2
        return result

    def _calc_in_sample_std_QA(self):
        for value in list(self._result_mapping.values()):
            flags = self.result_df[value + ' std'].apply(lambda x: 'ISV' if x > self._isv_thresholds[value] else '')
            self.result_df[value + ' QA'] = self._concat_qa_strings(self.result_df, value, flags)

    def _calc_bbv(self):
        for value in list(self._result_mapping.values()):
            # Only use unspiked samples for BBV calculation
            unspiked = self.result_df.loc[(slice(None), 0),:]
            # Group by bottle (Site Name, Sample Datetime)
            grouped = unspiked.set_index(['Site Name', 'Sample Datetime'])[[value + ' mean']].groupby(['Site Name', 'Sample Datetime'])
            # Calculate BBV flag
            flags = grouped.transform(lambda x: 'BBV' if (((x.max() - x.min())) > 2*self._isv_thresholds[value]) else '')
            flags = flags.fillna('')
            flags.index = unspiked.index
            unspiked[value + ' QA'] = self._concat_qa_strings(unspiked, value, flags[value + ' mean'])
            # Update the main result_df only for unspiked rows
            self.result_df.loc[unspiked.index, value + ' QA'] = unspiked[value + ' QA']

    def _check_pH(self):
        if 'Start pH after NaHCO3' in self.result_df.columns:
            flags = pd.DataFrame()
            flags.index = self.result_df.index
            flags['pH flag'] = ''
            flags[self.result_df['Start pH after NaHCO3']>self._pH_threshold] = 'PH'
            print(flags)
            for value in self._result_mapping.values():
                self.result_df[value + ' QA'] = self._concat_qa_strings(self.result_df, value, flags['pH flag'])



    def plot_drift_QA(self, result_name, **kwargs):
        drift_readings = self._raw_result_df[self._result_df['Cup Type'] == 'DRIF']
        drift_readings.plot(y=result_name, x='Date Time Stamp', **kwargs)

    def plot_QA(self, result_name, **kwargs):
        drift_readings = self._raw_result_df[self._raw_result_df['Cup Type'] == 'DRIF']
        baseline_readings = self._raw_result_df[self._raw_result_df['Cup Type'] == 'BASL']
        if 'ax' not in kwargs:
            fig, ax = plt.subplots()
            drift_readings.plot(y=result_name, x='Date Time Stamp', label = 'Drift', marker = '*', ax=ax, **kwargs)
            baseline_readings.plot(y=result_name, x='Date Time Stamp', label = 'Baseline', marker = 'o', ax=ax, **kwargs)
        else:
            drift_readings.plot(y=result_name, x='Date Time Stamp', label = 'Drift', marker = '*', **kwargs)
            baseline_readings.plot(y=result_name, x='Date Time Stamp', label = 'Baseline', marker = 'o',**kwargs)
