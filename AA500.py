import pandas as pd
import matplotlib.pyplot as plt

class AA500_Result:
    def __init__(self, result_path, samplelist_path, result_mapping={'Result 1':'Nitrate', 'Result 2':'Phosphate', 'Result 3':'Ammonium'}, isv_thresholds = {'Nitrate':.005, 'Phosphate':.004, 'Ammonium': .002}, bbv_thresholds = {'Nitrate':.1, 'Phosphate':.1, 'Ammonium': .1}):
        self._isv_thresholds = isv_thresholds
        self._bbv_thresholds = bbv_thresholds

        self._result_mapping = result_mapping
        self._raw_result_df = self._read_AA500_results(result_path)
        self._metadata = self._read_AA500_metadata(result_path)
        self._samplelist_df = self._read_samplelist_df(samplelist_path)
        self.result_df = self._merge_results_samplelist(result_mapping)
        self.result_df['AA500 Run Date'] = pd.to_datetime(self._metadata.loc['DATE', 'Value']+ ' ' + self._metadata.loc['TIME', 'Value'])
        self.result_df['AA500 Operator'] = self._metadata.loc['OPER', 'Value']
        self._calc_in_sample_std_QA()
        self._calc_bbv()

    def _read_AA500_metadata(self, path):
        metadata = pd.read_csv(path, nrows=7, names=['Param', 'Value'], usecols=[0,1]).set_index('Param')
        return metadata

    def _read_AA500_results(self, path):
        df = pd.read_csv(path, header=14, parse_dates=['Date Time Stamp']).rename(columns = self._result_mapping)
        return df

    def _read_samplelist_df(self, path):
        df = pd.read_excel(path)
        df['Sample Datetime'] = pd.to_datetime(df['Sample Collected'].astype(str) + ' ' + df['Time Sample Collected'].fillna('00:00:00').astype(str), errors='coerce')
        return df
    
    def _merge_results_samplelist(self, result_mapping):
        samples = self._raw_result_df[(self._raw_result_df['Cup Type'] == 'DSAMP')| (self._raw_result_df['Cup Type'] == '3SAMP') | (self._raw_result_df['Cup Type'] == 'SAMP')] 
        sample_results = samples.groupby('Cup Number')[list(result_mapping.values())].describe()
        sample_results= sample_results.loc[:, pd.IndexSlice[:, ['mean', 'std']]]

        err = sample_results.loc[:, pd.IndexSlice[:, 'std']].apply(lambda x: 2 * x)         ## add err (2* std)
        err.columns = pd.MultiIndex.from_tuples([(col[0], 'err') for col in err.columns])
        sample_results = pd.concat([sample_results, err], axis=1).sort_index(axis=1)

        sample_results.columns = sample_results.columns.map(' '.join)
        return self._samplelist_df.join(sample_results, on = 'Cup Number', how='right', validate='one_to_one').set_index(['Site Name', 'Sample Datetime']).sort_values(['Site Name', 'Sample Datetime'])
    
    def _concat_qa_strings(self, value, column2):
        result = pd.Series()
        if value + ' QA' in self.result_df.columns:
            column1  = self.result_df[value + ' QA']
            result = pd.Series([
                ','.join([a for a in [x, y] if a])
                for x, y in zip(column1, column2)
            ], index = self.result_df.index)
        else:
            result = column2
        return result

    def _calc_in_sample_std_QA(self):
        for value in list(self._result_mapping.values()):
            flags = self.result_df[value + ' std'].apply(lambda x: 'ISV' if x > self._isv_thresholds[value] else '')
            self.result_df[value + ' QA'] = self._concat_qa_strings(value, flags)

    def _calc_bbv(self):
        for value in list(self._result_mapping.values()):
            flags = self.result_df[[value+ ' mean']].reset_index().groupby(['Site Name', 'Sample Datetime']).transform(lambda x: 'BBV' if (((x.max()-x.min())/x.mean())> self._bbv_thresholds[value]) else '')
            flags = flags.fillna('')
            flags.index = self.result_df.index
            self.result_df[value + ' QA'] = self._concat_qa_strings(value, flags[value + ' mean'])

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
