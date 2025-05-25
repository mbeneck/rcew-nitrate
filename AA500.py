import pandas as pd
import matplotlib.pyplot as plt

class AA500_Result:
    def __init__(self, result_path, samplelist_path, result_mapping={'Result 1':'Nitrate', 'Result 2':'Phosphate', 'Result 3':'Ammonium'}, qa_thresholds = {'Nitrate':.005, 'Phosphate':.004, 'Ammonium': .002}):
        self._qa_thresholds = qa_thresholds
        self._result_mapping = result_mapping
        self._raw_result_df = self._read_AA500_results(result_path)
        self._samplelist_df = self._read_samplelist_df(samplelist_path)
        self.result_df = self._merge_results_samplelist(result_mapping)
        self._calc_in_sample_std_QA()
        

    def _read_AA500_results(self, path):
        df = pd.read_csv(path, header=14, parse_dates=['Date Time Stamp']).rename(columns = self._result_mapping)
        return df

    def _read_samplelist_df(self, path):
        df = pd.read_excel(path)
        df['Sample Datetime'] = pd.to_datetime(df['Sample Collected'].astype(str) + ' ' + df['Time Sample Collected'].astype(str), errors='coerce')
        return df
    
    def _merge_results_samplelist(self, result_mapping):
        samples = self._raw_result_df[(self._raw_result_df['Cup Type'] == 'DSAMP')| (self._raw_result_df['Cup Type'] == '3SAMP') | (self._raw_result_df['Cup Type'] == 'SAMP')] 
        sample_results = samples.groupby('Cup Number')[list(result_mapping.values())].describe()
        #sample_results.rename(result_mapping, axis=1, inplace=True)
        sample_results= sample_results.loc[:, pd.IndexSlice[:, ['mean', 'std']]]
        sample_results.columns = sample_results.columns.map(' '.join)
        return self._samplelist_df.join(sample_results, on = 'Cup Number', how='right', validate='one_to_one').set_index(['Site Name', 'Sample Datetime']).sort_values(['Site Name', 'Sample Datetime'])

    def _calc_in_sample_std_QA(self):
        for value in list(self._result_mapping.values()):
            print(self.result_df)
            self.result_df[value + ' QA'] = self.result_df[value + ' std'].apply(lambda x: 'ISV' if x > self._qa_thresholds[value] else '')


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

    def export_calibration(self, path):
        export_cal = pd.concat((self.intercepts/1000, self.slopes), axis=1)
        export_cal.T.to_csv(path, header=False, index=False, sep='\t')