def band_indices(sensorDict):
    import pandas as pd

    # OLCI
    s3rrs = sensorDict['s3']
    #index = ['SICF','SIPF','BAIR','PCI','FAI','BG','FT','b2','b3','ndci','pcS05','chlG05','G2b','G3b']
    #s3rrs_idx = pd.DataFrame(index=index)
    for i,k in enumerate(s3rrs.index):
        rrs = s3rrs.iloc[i,0:18]
        s3rrs.at[k,'SICF'] = rrs.loc[681.25] - rrs[665] - (rrs.loc[753.75] - rrs[665]) * (681.25 - 665) / (753.75 - 665)
        s3rrs.at[k,'SIPF'] = rrs[665] - rrs[620] - (rrs.loc[681.25] - rrs[620]) * (665 - 620) / (681.25 - 620)
        s3rrs.at[k,'BAIR'] = rrs.loc[708.75] - rrs.loc[681.25] - (rrs.loc[753.75] - rrs.loc[681.25]) * (708.75 - 681.25) / (753.75 - 681.25)
        s3rrs.at[k,'PCI'] = rrs[620] - rrs[560] - (rrs[665] - rrs[560]) * (620 - 560) / (665 - 560)
        s3rrs.at[k,'FAI'] = rrs.loc[753.75] - rrs.loc[681.25] - (rrs[885] - rrs.loc[681.25]) * (753.75 - 681.25) / (885 - 681.25)
        s3rrs.at[k,'BG'] = rrs[560] - rrs[510] - (rrs[620] - rrs[510]) * (560 - 510) / (620 - 510)
        s3rrs.at[k,'FT'] = rrs.loc[673.75] - rrs[665] - (rrs.loc[681.25] - rrs[665]) * (673.75 - 665) / (681.25 - 665)
        s3rrs.at[k,'b2'] = rrs.loc[708.75] / rrs[665]
        s3rrs.at[k,'b3'] = ((1/rrs[665] - 1/rrs.loc[708.75]) * rrs.loc[753.75])
        s3rrs.at[k,'ndci'] = (rrs.loc[708.75] - rrs[665]) / (rrs.loc[708.75] + rrs[665])

        # gons
        bb = (1.61 * rrs.loc[778.75]) / (0.082-0.6 * rrs.loc[778.75])
        chlaG05 = ((rrs.loc[708.75]/rrs[665])*(0.727+bb)-0.401-bb) / 0.015
        pcS05 = 170 * ((rrs.loc[708.75]/rrs[620])*(0.727+bb)-0.281-bb) - (0.51*chlaG05)
        s3rrs.at[k, 'pcS05'] = pcS05
        s3rrs.at[k, 'chlG05'] = chlaG05

        # Gilerson 2010
        G2b = (35.75 * (rrs.loc[708.75]/rrs[665]) - 19.3) ** 1.124
        G3b = (113.36 * ((1/rrs[665] - 1/rrs.loc[708.75]) * rrs.loc[753.75]) + 16.45) ** 1.124
        s3rrs.at[k, 'G2b'] = G2b
        s3rrs.at[k, 'G3b'] = G3b

    # MSI
    s2rrs = sensorDict['s2']
    #index = ['BAIR','FAI','BG','b2','b3','ndci']
    #s2rrs_idx = pd.DataFrame(index=index)
    for i,k in enumerate(s2rrs.index):
        rrs = s2rrs.iloc[i,0:9]
        s2rrs.at[k,'BAIR'] = (rrs[705] - rrs[665] - (rrs[740] - rrs[665]) * (705 - 665) / (740 - 665))
        s2rrs.at[k,'FAI'] = (rrs[783] - rrs[665] - (rrs[865] - rrs[665]) * (783 - 665) / (865 - 665))
        s2rrs.at[k,'BG'] = (rrs[560] - rrs[490] - (rrs[665] - rrs[490]) * (560 - 490) / (665 - 490))
        s2rrs.at[k,'b2'] = (rrs[705] / rrs[665])
        s2rrs.at[k,'b3'] = ((1/rrs[665] - 1/rrs[705]) * rrs[740])
        s2rrs.at[k,'ndci'] = (rrs[705] - rrs[665]) / (rrs[705] + rrs[665])

    # MERIS
    merisrrs = sensorDict['meris']
    # index = ['SICF','SIPF','BAIR','PCI','FAI','BG','b2','b3','ndci']
    # merisrrs_idx = pd.DataFrame(index=index)
    for i,k in enumerate(merisrrs.index):
        rrs = merisrrs.iloc[i,0:15]
        merisrrs.at[k,'SICF'] = rrs.loc[681.25] - rrs[665] - (rrs.loc[708.75] - rrs[665]) * (681.25 - 665) / (708.75 - 665)
        merisrrs.at[k,'SIPF'] = rrs[665] - rrs[620] - (rrs.loc[681.25] - rrs[620]) * (665 - 620) / (681.25 - 620)
        merisrrs.at[k,'BAIR'] = rrs.loc[708.75] - rrs.loc[681.25] - (rrs.loc[753.75] - rrs.loc[681.25]) * (708.75 - 681.25) / (753.75 - 681.25)
        merisrrs.at[k,'PCI'] = rrs[620] - rrs[560] - (rrs[665] - rrs[560]) * (620 - 560) / (665 - 560)
        merisrrs.at[k,'FAI'] = rrs.loc[753.75] - rrs.loc[681.25] - (rrs[885] - rrs.loc[681.25]) * (753.75 - 681.25) / (885 - 681.25)
        merisrrs.at[k,'BG'] = rrs[560] - rrs[510] - (rrs[620] - rrs[510]) * (560 - 510) / (620 - 510)
        merisrrs.at[k,'b2'] = rrs.loc[708.75] / rrs[665]
        merisrrs.at[k,'b3'] = ((1/rrs[665] - 1/rrs.loc[708.75]) * rrs.loc[753.75])
        merisrrs.at[k,'ndci'] = (rrs.loc[708.75] - rrs[665]) / (rrs.loc[708.75] + rrs[665])

    # MODIS
    modisrrs = sensorDict['modis']
    # index = ['SICF','SIPF','FAI','BG']
    # modisrrs_idx = pd.DataFrame(index=index)
    for idx in modisrrs.index:
        rrs = modisrrs.iloc[i,0:13]
        modisrrs.at[k,'SICF'] = rrs[678] - rrs[667] - (rrs[748] - rrs[667]) * (678 - 667) / (748 - 667)
        modisrrs.at[k,'SIPF'] = rrs[645] - rrs[555] - (rrs[678] - rrs[645]) * (645 - 555) / (678 - 555)
        modisrrs.at[k,'FAI'] = rrs[748] - rrs[678] - (rrs[869] - rrs[678]) * (748 - 678) / (869 - 678)
        modisrrs.at[k,'BG'] = rrs[555] - rrs[488] - (rrs[667] - rrs[488]) * (555 - 488) / (667 - 488)

    ##
    # info1 = s3rrs.iloc[:,-17:]
    # frames1 = [s3rrs_idx,info1]
    # frames2 = [s2rrs_idx,info1]
    # frames3 = [merisrrs_idx,info1]
    # frames4 = [modisrrs_idx,info1]
    # #
    # s3rrs_idx_final = pd.concat(frames1,axis=1)
    # s2rrs_idx_final = pd.concat(frames2,axis=1)
    # merisrrs_idx_final = pd.concat(frames3,axis=1)
    # modisrrs_idx_final = pd.concat(frames4,axis=1)
    #
    sensorDict['s3'] = s3rrs
    sensorDict['s2'] = s2rrs
    sensorDict['meris'] = merisrrs
    sensorDict['modis'] = modisrrs

    return sensorDict
