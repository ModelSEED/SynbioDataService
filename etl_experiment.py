import pandas as pd
from pandas import DataFrame
import numpy as np
import io
import re
from datetime import datetime


class AddPlateParameters:

    def __init__(self, protocol_id, operation_id, plate_type, plate_layout, measurement_type, lab, contact,
                 plate_id=None, plate_index=None, plate_transfer=None, plate_timestamp=None,
                 client=None):

        self.protocol_id = protocol_id
        self.operation_id = operation_id
        self.measurement_type = measurement_type
        self._lab = lab
        self._contact = contact
        self.timestamp = plate_timestamp
        self.index = plate_index
        self.transfer = plate_transfer
        self.id = plate_id
        self.type = plate_type
        self.layout = plate_layout
        if client:
            d_lab = pd.read_sql(f"SELECT * FROM lab WHERE id = {self._lab}", client).transpose().to_dict()
            d_people = pd.read_sql(f"SELECT * FROM people WHERE id = {self._contact}", client).transpose().to_dict()
            if len(d_lab) == 1:
                self._lab = d_lab[0]
            else:
                raise ValueError(f'Invalid Lab ID: {self._lab}')
            if len(d_people) == 1:
                self._contact = d_people[0]
            else:
                raise ValueError(f'Invalid Contact ID: {self._contact}')

    @property
    def lab_id(self):
        if type(self._lab) == dict:
            return self._lab['id']
        return self._lab

    @property
    def contact_id(self):
        if type(self._contact) == dict:
            return self._contact['id']
        return self._contact

    def _repr_html_(self):
        _str_lab = self._lab['name'] if type(self._lab) == dict else self._lab
        _str_contact = self._contact['email'] if type(self._contact) == dict else self._contact
        return """
        <table>
            <tr>
                <td><strong>ID</strong></td>
                <td>{address}</td>
            </tr><tr>
                <td><strong>Protocol</strong></td>
                <td>{protocol_id}</td>
            </tr><tr>
                <td><strong>Operation</strong></td>
                <td>{operation_id}</td>
            </tr><tr>
                <td><strong>measurement_type</strong></td>
                <td>{measurement_type}</td>
            </tr><tr>
                <td><strong>Lab</strong></td>
                <td>{contact} ({lab})</td>
            </tr><tr>
                <td><strong>plate</strong></td>
                <td>{plate_type} ({plate_layout})</td>
            </tr>
          </table>""".format(
            address="0x0%x" % id(self),
            protocol_id=self.protocol_id,
            operation_id=self.operation_id,
            measurement_type=self.measurement_type,
            contact=_str_contact,
            lab=_str_lab,
            plate_type=self.type,
            plate_layout=self.layout,
        )


class EtlExperiment:

    def __init__(self, engine, minio):
        self.engine = engine
        self.mio = minio

    def get_plate_reader_filenames(self, minio_bucket, minio_path_to_plate_reader_files, regexpattern):
        # Return list of file names that match the expected file name pattern.
        objects = self.mio.list_objects(minio_bucket, prefix=minio_path_to_plate_reader_files)

        filepaths = [obj.object_name for obj in objects]
        filenames = [fn.split('/')[1] for fn in filepaths]
        plate_reader_files = [fn for fn in filenames if re.match(regexpattern, fn)]

        return plate_reader_files

    def load_layout_file(self, minio_bucket_name, layout_file_path):
        # Read plate layout file.
        response = None
        df = None
        try:
            response = self.mio.get_object(minio_bucket_name, layout_file_path)
            csv_data = response.data
            df = pd.read_csv(io.BytesIO(csv_data), header=None)
        except Exception as e:
            print(f"Error: {e}")
        finally:
            if response:
                response.close()
                response.release_conn()
            return df

    def load_layout(self, minio_bucket_name, path_to_layout_files):
        transfer_layout = self.load_layout_file(minio_bucket_name, path_to_layout_files + 'transfer_layout.csv')
        rep_layout = self.load_layout_file(minio_bucket_name, path_to_layout_files + 'replicate_layout.csv')
        strain_layout = self.load_layout_file(minio_bucket_name, path_to_layout_files + 'strain_layout.csv')
        gc_layout = self.load_layout_file(minio_bucket_name, path_to_layout_files + 'growth_condition_layout.csv')
        return transfer_layout, rep_layout, strain_layout, gc_layout

    @staticmethod
    def get_transfers(batch_dict, batch, plate, transfer):
        # Which transfer since beginning of experiment?
        return batch_dict[batch] + (plate - 1) * 3 + transfer

    @staticmethod
    def get_plates(batch_dict, batch, plate):
        # How many plates in total?
        return batch_dict[batch] + plate

    @staticmethod
    def get_parent_plate(plate, column):
        # Which plate was the parent sample on?
        if column in (1, 4, 7):
            parent_plate = plate - 1
        else:
            parent_plate = plate
        return parent_plate

    @staticmethod
    def get_parent_col(col):
        # Which column was the parent sample in?
        cols = (1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 10, 11)
        parent_cols = (3, 1, 2, 6, 4, 5, 9, 7, 8, 0, 0, 0)
        parent_col_dict = dict(zip(cols, parent_cols))
        return parent_col_dict[col]

    @staticmethod
    def get_well_name(row, col):
        # Return well name based on plate rows/cols.
        well_name = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'][row] + str(col + 1)
        return well_name

    def etl_plate(self, plate_data: DataFrame, plate_parameters: AddPlateParameters, start_date,
                  exp_id, exp_index, minio_bucket_name='synbio'):
        """
        Extracts and transforms plate reader data from MinIO storage,
        and stores it in MySQL database.
    
        Args:
    
        Returns:
        """

        transfer_layout, rep_layout, strain_layout, gc_layout = self.load_layout(minio_bucket_name,
                                                                                 plate_parameters.layout)

        #plate_reader_filenames = self.get_plate_reader_filenames(minio_bucket_name,
        #                                                         path_to_plate_reader_files,
        #                                                         fname_pattern)
        
        # Initialize data df
        data = pd.DataFrame()

        # Parse info contained in plate reader file name
        data_row = {}
        data_row['file_ID'] = plate_parameters.id
        data_row['timestamp'] = plate_parameters.timestamp
        data_row['datetime'] = datetime.fromtimestamp(data_row['timestamp']).isoformat()
        data_row['plate_index'] = plate_parameters.index
        # t_transfer indicates which plate cols were most recently innoculated.
        data_row['t_transfer'] = plate_parameters.transfer

        # Read plate reader files into dataframe
        for row in range(8):
            for col in range(12):
                data_row['row'] = row
                data_row['column'] = col
                data_row['OD'] = plate_data.iloc[row, col]
                data = pd.concat([data, pd.Series(data_row).to_frame().T])
        
        data.reset_index(inplace=True, drop=True)
        
        # Translate row and column numbers to well names
        data['well'] = data.apply(
            lambda x: self.get_well_name(x['row'], x['column']), axis=1
            )
        # Translate timestamp into isoformat time
        data['datetime'] = data['timestamp'].apply(
            lambda x: datetime.fromtimestamp(x).isoformat()
            )
        
        # Appending metadata
        # data['filename'] = path_to_plate_reader_files
        data['measurement_type'] = plate_parameters.measurement_type
        data['experiment'] = exp_id
        data['plate_type'] = plate_parameters.type
        data['start_date'] = start_date
        data['exp_index'] = exp_index
        data['operation_id'] = plate_parameters.operation_id
        data['layout_filename'] = plate_parameters.layout
        
        # Given location on plate (row, col) and layout files, get the strain, 
        # growth condition, replicate number, and transfer_l for each well.
        # (transfer_l indicates transfer based on plate location.)
        data['strain'] = data.apply(
            lambda x: strain_layout.iloc[x['row'], x['column']], axis=1
            ).astype("Int64")
        data['replicate'] = data.apply(
            lambda x: rep_layout.iloc[x['row'], x['column']], axis=1
            ).astype("Int64")
        data['gc'] = data.apply(
            lambda x: gc_layout.iloc[x['row'], x['column']], axis=1
            ).astype("Int64")
        data['l_transfer'] = data.apply(
            lambda x: transfer_layout.iloc[x['row'], x['column']], axis=1
            ).astype("Int64")
        
        # Determine innoculation status
        # Explanation: All wells on a plate are read at each timepoint, but only some
        # of the wells will be inoculated at any given timepoint. Only wells that have
        # media in them (i.e., growth condition not NA) are considered *samples*.
        # Samples that do not have a strain designation are negative controls.
        # Depending on the location on the plate, samples can have 11, 22, or 33
        # timepoints. The calculation below determines which transfer to assign a
        # given OD reading to.
        
        # Calculate actual transfer number for a given reading
        data['l-t_transfer'] = data.apply(
            lambda x: x['l_transfer']-x['t_transfer'], axis=1
            )
        data['transfer'] = np.where(
            (data['l-t_transfer']<=0), data['l_transfer'], 1
            )
        
        # Compute total number of plates, transfers for this experiment
        # There are 3 transfers per plate, x plates per batch, and y batches per
        # experiment. A batch is defined as a continuous run of the robot without
        # any interruption. All measurements in a batch will have the same file_ID.
        # Each batch starts with plate=1, transfer=1, timepoint=1. To accurately
        # calculate the cummulative number of transfers since the beginning of the
        # experiment for each data point, we need to calculate the number of transfers
        # for the individual batches and their cumulative number (in their correct 
        # order).
        
        batches = (
            data.groupby('file_ID')['datetime']
            .agg("min")
            .reset_index()
            .sort_values('datetime')['file_ID']
            .to_list()
        )
        batch_n_transfers = [0]  # transfers for each of the batches.
        # Start with 0, because the first batch has no prior transfers.
        batch_n_plates = [0]  # number of plates for each of the batches.
        # Start with 0, because first batch has no prior plates.
        
        for batch in batches:
            # Loop through each batch and calculate the number of transfers
            batch_data = data.loc[
                (data['file_ID'] == batch) & (~pd.isna(data['transfer']))
            ]
            # From last measurement in batch, total number of plates and transfers in
            # this batch
            max_plate = batch_data.sort_values(['plate_index', 'transfer']).iloc[-1]
            batch_n_plates.append(max_plate['plate_index'])
            batch_transfers = (max_plate['plate_index']-1)*3 + (max_plate['transfer'])
            batch_n_transfers.append(batch_transfers)
        
        # cumulative transfers (# of transfers since the start of the experiment)
        batch_cumsum_t = np.cumsum(np.array(batch_n_transfers[:-1])) 
        batch_dict_t = dict(zip(batches, batch_cumsum_t))
        # cumulative plates (# of plates since the start of the experiment)
        batch_cumsum_p = np.cumsum(np.array(batch_n_plates[:-1])) 
        batch_dict_p = dict(zip(batches, batch_cumsum_p))
        
        data['cum_transfer'] = np.where(
            data['transfer']==0, np.array([0]*len(data)), 
            data.apply(
                lambda x: self.get_transfers(
                    batch_dict_t, x['file_ID'], x['plate_index'], x['transfer']
                    ), axis=1)
                    ) 
        data['cum_plate'] = data.apply(
            lambda x: self.get_plates(
                batch_dict_p, x['file_ID'], x['plate_index']
                ), axis=1
                )
        
        # data['strain'] is NA (i.e., negative control) if not yet inoculated
        data['strain'] = np.where(
            (data['l-t_transfer']>0), 
            np.array([pd.NA]*len(data)), 
            data['strain']
            )
        data['strain'] = data['strain'].astype("Int64")
        data['replicate'] = np.where(
            pd.isna(data['strain']), 
            np.array([pd.NA]*len(data)), 
            data['replicate']
            )
        data['replicate'] = data['replicate'].astype("Int64")
        
        # Calculate a background value for each plate reader measurement
        # (based on the wells that only contain media)
        data['background'] = pd.NA
        data['background'] = data.groupby(
            ['experiment', 'plate_index', 'timestamp']
            )['OD'].transform(
            lambda x: x[data.loc[x.index, 'strain'].isna()].mean()
            )
        
        # Compute innoculation timestamp based on the oldest timestamp
        data['innoculation_timestamp'] = pd.NA
        data.loc[
            (~(pd.isna(data['cum_transfer']))) & (~(pd.isna(data['strain']))),
            'innoculation_timestamp'
            ] = data.groupby(
                ['cum_plate', 'cum_transfer']
                )['datetime'].transform("min")
        
        # Assign parent samples
        # Only innoculated samples (i.e., not neg. controls) can have parent samples.
        # Only samples after passage 1 have parent samples (passage 1 parents will 
        # have to be manually assigned)
        data['parent_plate'] = pd.NA
        data.loc[
            ~pd.isna(data['strain']) & (data['cum_transfer'] != 1),
            'parent_plate'] = data.apply(
                lambda x: self.get_parent_plate(
                    x['cum_plate'], x['column']
                    ), axis=1)
        data['parent_well'] = pd.NA
        data.loc[
            ~pd.isna(data['strain']) & (data['cum_transfer'] != 1),
            'parent_well'] = data.apply(
                lambda x: self.get_well_name(x['row'], self.get_parent_col(x['column'])), axis=1)
        
        # Assign plate, well, and sample names
        data = data.convert_dtypes()
        data['plate_name'] = (
            'E:' + data['experiment'] + '.P:' + data['cum_plate'].astype(str)
        )
        data['well_name'] = data['plate_name'] + '.W:' + data['well'].astype(str)
        data['sample_name']= pd.NA
        data.loc[(~(pd.isna(data['gc']))), 'sample_name'] = (
            data['well_name'] + '.S:' + data['strain'].astype(str) + '.C:' +
            data['gc'].astype(str) + '.R:' + data['replicate'].astype(str) +
            '.T:' + data['cum_transfer'].astype(str)
        )
        
        # get parent_id
        data['parent_id'] = pd.NA
        data.loc[
            (~(pd.isna(data['strain']))) &
            (~(pd.isna(data['gc']))) &
            (data['cum_transfer'] != 1),
            'parent_id'] = (
                'E:' + data['experiment'] + '.P:' + data['parent_plate'].astype(str) + 
                '.W:' + data['parent_well'].astype(str) + '.S:' + data['strain'].astype(str) + 
                '.C:' + data['gc'].astype(str) + '.R:' + data['replicate'].astype(str) + 
                '.T:' + (data['cum_transfer']-1).astype(str)
            )
        # data = data.convert_dtypes()
        
        
        # Reformat data for database upload
        # 1) Plates
        plates = data[
            ['plate_name', 'experiment', 'plate_type', 'plate_index', 'layout_filename']
            ].drop_duplicates().reset_index(drop=True)
        plates.rename(
            columns={'plate_name': 'id', 'experiment': 'experiment_id'}, inplace=True
            )
        plates.to_sql('plate', self.engine, index=False, if_exists='append')
        
        # 2) Samples and associated measurements
        # Each sample has a measurement of type 'growth' to which all od_measurements map.
        sample_meas = data.loc[
            ~pd.isna(data['sample_name']) & (~(pd.isna(data['gc'])))
            ].drop_duplicates(
                (['sample_name', 'experiment', 'plate_name', 'well', 'cum_transfer',
                  'gc', 'strain', 'innoculation_timestamp', 'replicate', 'parent_id'])
                         ).sort_values('innoculation_timestamp')
        samples = sample_meas[
            (['sample_name', 'experiment', 'plate_name', 'well', 'cum_transfer','gc', 
              'strain', 'innoculation_timestamp', 'replicate', 'parent_id'])
              ].copy()
        samples.rename(
            columns={
                'sample_name': 'name', 'experiment': 'experiment_id', 
                'plate_name': 'plate', 'cum_transfer': 'passage',
                'gc': 'growth_condition_id', 'strain': 'strain_id',
                'parent_id': 'parent_sample_name'
                }, inplace=True)
        samples.to_sql('sample', self.engine, index=False, if_exists='append')
        measurements = sample_meas[
            ['sample_name', 'operation_id', 'measurement_type', 'filename']
            ].copy()
        measurements.rename(
            columns={'sample_name': 'sample_id', 'measurement_type': 'type'},
            inplace=True
            )
        measurements.to_sql('measurement', self.engine, index=False, if_exists='append')
        
        # 3) OD_measurements
        # The index of the db table `measurements` autoincrements, so the measurement
        # ids to which the od_measurements of the current dataset belong are not known
        # ahead of time. Therefore, need to download all measurements whith sample ids
        # from this dataset to get their ids, then merge with OD readings.
        sample_names = tuple(sample_meas['sample_name'])
        meas_from_db = pd.read_sql(
            f"SELECT `id`, `sample_id` FROM `measurement` WHERE `sample_id` IN {sample_names};",
            self.engine
            ).rename(columns={'id': 'measurement_id'})
        od_meas = meas_from_db.merge(
            data, left_on='sample_id', right_on='sample_name', how='inner'
            )[['plate_name', 'measurement_id', 'datetime', 'OD', 'background']].rename(
                columns={'OD': 'od'}
                )
        od_meas.drop('plate_name', axis=1).to_sql('od_measurement', self.engine, index=False, if_exists='append')
    
        
        # 4) Generate growth curves with AMiGA code and store in db   
        # Calculate background value, calculate timepoints in hours.
        
        # First drop duplicate measurement.
        # (There are always two timepoints that are the same in the ALE1a dataset
        # because measurement taken before and after transfer received same
        # timepoint.)
        od_meas = od_meas.drop_duplicates(['measurement_id', 'datetime']).copy()
        od_meas['od_sub_background'] = od_meas['od'] - od_meas['background']
        
        # For each sample (measurement_id), calculate timepoints from datetime in hours
        od_meas['datetime'] = od_meas['datetime'].apply(lambda x: datetime.fromisoformat(x))
        baseline_timepoints = od_meas.groupby('measurement_id')['datetime'].transform("min")
        od_meas['timepoint'] = od_meas['datetime'] - baseline_timepoints
        od_meas['timepoint'] = od_meas['timepoint'].apply(lambda x: x.total_seconds()/3600)
        
        # Collection data frame
        metrics_all = pd.DataFrame()
        
        # For every plate_name, create a GrowthPlate
        
        value = 'od' #'OD_sub_background' # or
        
        for plate_name, plate_data in od_meas.groupby('plate_name'):
        
            pivot_data = pd.pivot(plate_data, 
                                  index = 'timepoint',
                                  values = value,
                                  columns = 'measurement_id'
                                 ).reset_index(
                                 ).rename(columns=({'timepoint':'time'})
                                 ).sort_values('time'
                                 ).astype(float)
            mapping = pd.DataFrame(data = {
                'm_id':plate_data['measurement_id'].unique(),
                'plate_ID':[plate_name]*len(plate_data['measurement_id'].unique())
            }
                                  ).set_index('m_id')
        
            plate = GrowthPlate(data=pivot_data,key=mapping)
            plate.computeBasicSummary()
            plate.raiseData()
            plate.logData()
            plate.subtractBaseline(True,poly=False)
            
            metrics = plate.key[['OD_Max']]
            metrics.insert(0, 'growth_rate', [pd.NA]*len(metrics))
            metrics.insert(0, 'doubling_time', [pd.NA]*len(metrics))
            metrics.insert(0, 'lag_time', [pd.NA]*len(metrics))
        
            for sample in plate.data.columns:
                # fit curve
                thisSample = pd.concat([plate.time, plate.data[sample]], axis=1).dropna()
                thisSample.columns = ['Time', 'OD']
                gm = GrowthModel(thisSample,ARD=True)
                gm.fit()
                curve = gm.run()
                metrics.loc[sample, 'growth_rate'] = curve.gr
                metrics.loc[sample, 'doubling_time'] = curve.td * 60 # convert to minutes
                metrics.loc[sample, 'lag_time'] = curve.lagC * 60 # convert to minutes
        
            metrics_all = pd.concat([metrics_all, metrics])
        
        
        # Massage data so it will fit into database
        # i.e., replace infinity values, replace very large values with NA,
        # and round to 3 dec places, 
        
        metrics_all = metrics_all.reset_index().rename(columns={
            'm_id':'measurement_id',
            'OD_Max': 'max_od'}
                                                      ).copy()
        metrics_all['lag_time'] = metrics_all['lag_time'].replace([np.inf, -np.inf], 10000)
        metrics_all['doubling_time'] = metrics_all['doubling_time'].replace([np.inf, -np.inf], 10000)
        metrics_all['max_od'] = metrics_all['max_od'].apply(lambda x: round(x, 3))
        metrics_all['growth_rate'] = metrics_all['growth_rate'].apply(lambda x: round(x, 3))
        
        crazy_lag_time = (metrics_all['lag_time'] >= 10000) | (metrics_all['lag_time'] < 0)
        metrics_all['lag_time'] = np.where(
            crazy_lag_time,
            [pd.NA]*len(metrics_all),
            metrics_all['lag_time'].apply(round)
        )
        
        crazy_doubling_time = (metrics_all['doubling_time'] >= 10000) | (metrics_all['doubling_time'] < 0)
        metrics_all['doubling_time'] = np.where(
            crazy_doubling_time,
            [pd.NA]*len(metrics_all),
            metrics_all['doubling_time'].apply(round)
        )
        
        od_max_near_0 = metrics_all['max_od'] < 0.0001
        metrics_all['max_od'] = np.where(
            od_max_near_0,
            [pd.NA]*len(metrics_all),
            metrics_all['max_od']
        )
        growth_rate_near_0 = metrics_all['growth_rate'] < 0.0001
        metrics_all['growth_rate'] = np.where(
            growth_rate_near_0,
            [pd.NA]*len(metrics_all),
            metrics_all['growth_rate']
        )
        
        metrics_all.to_sql('growth_measurement', self.engine, index=False, if_exists='append')

