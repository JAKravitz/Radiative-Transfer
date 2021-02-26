def sample_info(flPath,noFlPath,sensorDict):
    import pandas as pd
    import numpy as np
    import os

    wv = np.arange(400, 901, 1)
    for k in sensorDict['meris'].index:

        info = k.split('_')
        chl = float(info[0])
        admix = float(info[1])
        dinoD = float(info[2])
        cnap = float(info[3])
        cdom = float(info[4])
        theta = int(info[5])
        if len(info) == 7:
            cy = info[6]
        else:
            cy = 'nan'
        fqy = np.nan

        # fl amplitude calc
        noFl = '_'.join(info)
        noFlname = noFl+'.csv'

        flData = pd.read_csv(os.path.join(flPath, k+'.csv'), index_col=0)
        flData = np.interp(wv, flData.index, flData.Rrs)

        nofldata = pd.read_csv(os.path.join(noFlPath,noFlname),index_col=0)
        nofldata = np.interp(wv, nofldata.index, nofldata.Rrs)

        # fluorescence amplitude
        fl_amp = ((flData[285] - nofldata[285]) * .993) / .001

        # other sample info
        #nn_name1 = '_'.join(info[:-1])
        nn_name = 'NN_' + k + '.csv'
        nn_info = pd.read_csv('~/.wine/drive_c/HE52_data/NN/cyano/{0}'.format(nn_name),index_col=0).T
        #nn_info = pd.read_csv('/Volumes/500gb/hy_data/NN_files/{0}'.format(nn_name), index_col=0).T

        # info for meris
        merisrrs = sensorDict['meris']
        merisrrs.at[k,'chl'] = chl
        merisrrs.at[k,'admix'] = admix
        merisrrs.at[k,'dinoD'] = dinoD
        merisrrs.at[k,'cnap'] = cnap
        merisrrs.at[k,'cdom'] = cdom
        merisrrs.at[k,'FQY'] = fqy
        merisrrs.at[k,'cy'] = cy
        merisrrs.at[k,'fl_amp'] = fl_amp
        merisrrs.at[k,'pcS'] = nn_info['PC_simis'].values
        merisrrs.at[k,'pcL'] = nn_info['PC_li'].values
        merisrrs.at[k,'aphy440'] = nn_info['aphy440'].values
        merisrrs.at[k, 'ag440'] = nn_info['ag440'].values
        merisrrs.at[k, 'anap440'] = nn_info['anap440'].values
        merisrrs.at[k, 'bphy440'] = nn_info['bphy440'].values
        merisrrs.at[k, 'bnap440'] = nn_info['bnap440'].values
        merisrrs.at[k, 'bbphy440'] = nn_info['bbphy440'].values
        merisrrs.at[k, 'bbnap440'] = nn_info['bbnap440'].values

        # modis
        modisrrs = sensorDict['modis']
        modisrrs.at[k,'chl'] = chl
        modisrrs.at[k,'admix'] = admix
        modisrrs.at[k,'dinoD'] = dinoD
        modisrrs.at[k,'cnap'] = cnap
        modisrrs.at[k,'cdom'] = cdom
        modisrrs.at[k,'FQY'] = fqy
        modisrrs.at[k,'cy'] = cy
        modisrrs.at[k,'fl_amp'] = fl_amp
        modisrrs.at[k,'pcS'] = nn_info['PC_simis'].values
        modisrrs.at[k,'pcL'] = nn_info['PC_li'].values
        modisrrs.at[k,'aphy440'] = nn_info['aphy440'].values
        modisrrs.at[k, 'ag440'] = nn_info['ag440'].values
        modisrrs.at[k, 'anap440'] = nn_info['anap440'].values
        modisrrs.at[k, 'bphy440'] = nn_info['bphy440'].values
        modisrrs.at[k, 'bnap440'] = nn_info['bnap440'].values
        modisrrs.at[k, 'bbphy440'] = nn_info['bbphy440'].values
        modisrrs.at[k, 'bbnap440'] = nn_info['bbnap440'].values


        # l8
        l8rrs = sensorDict['l8']
        l8rrs.at[k,'chl'] = chl
        l8rrs.at[k,'admix'] = admix
        l8rrs.at[k,'dinoD'] = dinoD
        l8rrs.at[k,'cnap'] = cnap
        l8rrs.at[k,'cdom'] = cdom
        l8rrs.at[k,'FQY'] = fqy
        l8rrs.at[k,'cy'] = cy
        l8rrs.at[k,'fl_amp'] = fl_amp
        l8rrs.at[k,'pcS'] = nn_info['PC_simis'].values
        l8rrs.at[k,'pcL'] = nn_info['PC_li'].values
        l8rrs.at[k,'aphy440'] = nn_info['aphy440'].values
        l8rrs.at[k, 'ag440'] = nn_info['ag440'].values
        l8rrs.at[k, 'anap440'] = nn_info['anap440'].values
        l8rrs.at[k, 'bphy440'] = nn_info['bphy440'].values
        l8rrs.at[k, 'bnap440'] = nn_info['bnap440'].values
        l8rrs.at[k, 'bbphy440'] = nn_info['bbphy440'].values
        l8rrs.at[k, 'bbnap440'] = nn_info['bbnap440'].values

        # hico
        hicorrs = sensorDict['hico']
        hicorrs.at[k,'chl'] = chl
        hicorrs.at[k,'admix'] = admix
        hicorrs.at[k,'dinoD'] = dinoD
        hicorrs.at[k,'cnap'] = cnap
        hicorrs.at[k,'cdom'] = cdom
        hicorrs.at[k,'FQY'] = fqy
        hicorrs.at[k,'cy'] = cy
        hicorrs.at[k,'fl_amp'] = fl_amp
        hicorrs.at[k,'pcS'] = nn_info['PC_simis'].values
        hicorrs.at[k,'pcL'] = nn_info['PC_li'].values
        hicorrs.at[k,'aphy440'] = nn_info['aphy440'].values
        hicorrs.at[k, 'ag440'] = nn_info['ag440'].values
        hicorrs.at[k, 'anap440'] = nn_info['anap440'].values
        hicorrs.at[k, 'bphy440'] = nn_info['bphy440'].values
        hicorrs.at[k, 'bnap440'] = nn_info['bnap440'].values
        hicorrs.at[k, 'bbphy440'] = nn_info['bbphy440'].values
        hicorrs.at[k, 'bbnap440'] = nn_info['bbnap440'].values


        # olci
        s3rrs = sensorDict['s3']
        s3rrs.at[k,'chl'] = chl
        s3rrs.at[k,'admix'] = admix
        s3rrs.at[k,'dinoD'] = dinoD
        s3rrs.at[k,'cnap'] = cnap
        s3rrs.at[k,'cdom'] = cdom
        s3rrs.at[k,'FQY'] = fqy
        s3rrs.at[k,'cy'] = cy
        s3rrs.at[k,'fl_amp'] = fl_amp
        s3rrs.at[k,'pcS'] = nn_info['PC_simis'].values
        s3rrs.at[k,'pcL'] = nn_info['PC_li'].values
        s3rrs.at[k,'aphy440'] = nn_info['aphy440'].values
        s3rrs.at[k, 'ag440'] = nn_info['ag440'].values
        s3rrs.at[k, 'anap440'] = nn_info['anap440'].values
        s3rrs.at[k, 'bphy440'] = nn_info['bphy440'].values
        s3rrs.at[k, 'bnap440'] = nn_info['bnap440'].values
        s3rrs.at[k, 'bbphy440'] = nn_info['bbphy440'].values
        s3rrs.at[k, 'bbnap440'] = nn_info['bbnap440'].values


        # s2
        s2rrs = sensorDict['s2']
        s2rrs.at[k,'chl'] = chl
        s2rrs.at[k,'admix'] = admix
        s2rrs.at[k,'dinoD'] = dinoD
        s2rrs.at[k,'cnap'] = cnap
        s2rrs.at[k,'cdom'] = cdom
        s2rrs.at[k,'FQY'] = fqy
        s2rrs.at[k,'cy'] = cy
        s2rrs.at[k,'fl_amp'] = fl_amp
        s2rrs.at[k,'pcS'] = nn_info['PC_simis'].values
        s2rrs.at[k,'pcL'] = nn_info['PC_li'].values
        s2rrs.at[k,'aphy440'] = nn_info['aphy440'].values
        s2rrs.at[k, 'ag440'] = nn_info['ag440'].values
        s2rrs.at[k, 'anap440'] = nn_info['anap440'].values
        s2rrs.at[k, 'bphy440'] = nn_info['bphy440'].values
        s2rrs.at[k, 'bnap440'] = nn_info['bnap440'].values
        s2rrs.at[k, 'bbphy440'] = nn_info['bbphy440'].values
        s2rrs.at[k, 'bbnap440'] = nn_info['bbnap440'].values

        # hyper
        hyperrs = sensorDict['hyper']
        hyperrs.at[k,'chl'] = chl
        hyperrs.at[k,'admix'] = admix
        hyperrs.at[k,'dinoD'] = dinoD
        hyperrs.at[k,'cnap'] = cnap
        hyperrs.at[k,'cdom'] = cdom
        hyperrs.at[k,'FQY'] = fqy
        hyperrs.at[k,'cy'] = cy
        hyperrs.at[k,'fl_amp'] = fl_amp
        hyperrs.at[k,'pcS'] = nn_info['PC_simis'].values
        hyperrs.at[k,'pcL'] = nn_info['PC_li'].values
        hyperrs.at[k,'aphy440'] = nn_info['aphy440'].values
        hyperrs.at[k, 'ag440'] = nn_info['ag440'].values
        hyperrs.at[k, 'anap440'] = nn_info['anap440'].values
        hyperrs.at[k, 'bphy440'] = nn_info['bphy440'].values
        hyperrs.at[k, 'bnap440'] = nn_info['bnap440'].values
        hyperrs.at[k, 'bbphy440'] = nn_info['bbphy440'].values
        hyperrs.at[k, 'bbnap440'] = nn_info['bbnap440'].values

    sensorDict2 = {'l8': l8rrs, 'modis': modisrrs, 'meris': merisrrs, 's3': s3rrs, 's2': s2rrs, 'hico': hicorrs,
                  'hyper': hyperrs}

    return sensorDict2

