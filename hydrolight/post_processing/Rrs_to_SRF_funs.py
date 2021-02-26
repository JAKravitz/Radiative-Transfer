def resample_OLCI(rrsData,srfpath):
    import numpy as np
    import pandas as pd
    # srfs
    srfs = pd.read_csv(srfpath,index_col=0)
    srfs = srfs.replace(-999, 0)

    # integrate
    srfInt = np.trapz(srfs, x=srfs.index, axis=0)

    # extract band signal
    wl = np.array([400,412.5,442.5,490,510,560,620,665,673.75,681.25,708.75,753.75,761.25,764.375,767.75,778.75,865,885])
    chanRrs = np.zeros((len(wl),rrsData.shape[1]))

    for ichan in np.arange(len(wl)):

        chanFilter = np.tile(srfs.iloc[:,ichan], [rrsData.shape[1],1])
        chanFilter = chanFilter.T
        chanLProduct = chanFilter * rrsData
        chanLProduct.index = srfs.index
        chanInt = np.trapz(chanLProduct, x=chanLProduct.index, axis=0) #/ srfInt[ichan]
        chanRrs[ichan,:] = [x / srfInt[ichan] for x in chanInt]

    chanRrs = pd.DataFrame(chanRrs, index=wl)
    return chanRrs

def resample_s2(rrsData,srfpath):
    import numpy as np
    import pandas as pd
    # srfs
    srfs = pd.read_csv(srfpath,index_col=0)
    #srfs = srfs.replace(-999, 0)

    # integrate
    srfInt = np.trapz(srfs, x=srfs.index, axis=0)

    # extract band signal
    wl = np.array([443,490,560,665,705,740,783,842,865])
    chanRrs = np.zeros((len(wl),rrsData.shape[1]))

    for ichan in np.arange(len(wl)):

        chanFilter = np.tile(srfs.iloc[:,ichan], [rrsData.shape[1],1])
        chanFilter = chanFilter.T
        chanLProduct = chanFilter * rrsData
        chanLProduct.index = srfs.index
        chanInt = np.trapz(chanLProduct, x=chanLProduct.index, axis=0) #/ srfInt[ichan]
        chanRrs[ichan,:] = [x / srfInt[ichan] for x in chanInt]

    chanRrs = pd.DataFrame(chanRrs, index=wl)
    return chanRrs

def resample_l8(rrsData,srfpath):
    import numpy as np
    import pandas as pd
    # srfs
    srfs = pd.read_csv(srfpath,index_col=0)
    #srfs = srfs.replace(-999, 0)

    # integrate
    srfInt = np.trapz(srfs, x=srfs.index, axis=0)

    # extract band signal
    wl = np.array([443,482,561,665,865])
    chanRrs = np.zeros((len(wl),rrsData.shape[1]))

    for ichan in np.arange(len(wl)):

        chanFilter = np.tile(srfs.iloc[:,ichan], [rrsData.shape[1],1])
        chanFilter = chanFilter.T
        chanLProduct = chanFilter * rrsData
        chanLProduct.index = srfs.index
        chanInt = np.trapz(chanLProduct, x=chanLProduct.index, axis=0) #/ srfInt[ichan]
        chanRrs[ichan,:] = [x / srfInt[ichan] for x in chanInt]

    chanRrs = pd.DataFrame(chanRrs, index=wl)
    return chanRrs

def resample_modis(rrsData,srfpath):
    import numpy as np
    import pandas as pd
    # srfs
    srfs = pd.read_csv(srfpath,index_col=0)
    srfs = srfs.dropna()

    # integrate
    srfInt = np.trapz(srfs, x=srfs.index, axis=0)

    # extract band signal
    wl = np.array([412,443,469,488,531,551,555,645,667,678,748,859,869])
    chanRrs = np.zeros((len(wl),rrsData.shape[1]))

    for ichan in np.arange(len(wl)):

        chanFilter = np.tile(srfs.iloc[:,ichan], [rrsData.shape[1],1])
        chanFilter = chanFilter.T
        chanLProduct = chanFilter * rrsData
        chanLProduct.index = srfs.index
        chanInt = np.trapz(chanLProduct, x=chanLProduct.index, axis=0) #/ srfInt[ichan]
        chanRrs[ichan,:] = [x / srfInt[ichan] for x in chanInt]

    chanRrs = pd.DataFrame(chanRrs, index=wl)
    return chanRrs

def resample_meris(rrsData,srfpath):
    import numpy as np
    import pandas as pd
    # srfs
    srfs = pd.read_csv(srfpath,index_col=0)
    srfs = srfs.dropna()

    # integrate
    srfInt = np.trapz(srfs, x=srfs.index, axis=0)

    # extract band signal
    wl = np.array([412.5,442.5,490,510,560,620,665,681.25,708.75,753.75,760.625,778.75,865,885,900])
    chanRrs = np.zeros((len(wl),rrsData.shape[1]))

    for ichan in np.arange(len(wl)):

        chanFilter = np.tile(srfs.iloc[:,ichan], [rrsData.shape[1],1])
        chanFilter = chanFilter.T
        chanLProduct = chanFilter * rrsData
        chanLProduct.index = srfs.index
        chanInt = np.trapz(chanLProduct, x=chanLProduct.index, axis=0) #/ srfInt[ichan]
        chanRrs[ichan,:] = [x / srfInt[ichan] for x in chanInt]

    chanRrs = pd.DataFrame(chanRrs, index=wl)
    return chanRrs