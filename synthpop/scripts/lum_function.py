import numpy as np
import matplotlib.pyplot as plt
import pandas
import pdb


def apply_extinction(df, R_V=3.1):
    # filter_keys = ['unknown','J','Z087','H','W149','W169']
    filter_keys = ['Z087', 'Y106', 'J129', 'H158', 'W146', 'F184']
    eff_wavelengths = {'Z087': 0.863390, 'Y106': 1.045809, 'J129': 1.274700,
        'H158': 1.558247, 'W146': 1.301380, 'F184': 1.830067}
    # eff_wavelengths = {'unknown':0.863390,'J':1.228538,'Z087':0.863390,'H':1.638610,'W149':1.301380,'W169':1.301380}

    # calc A_lambda/A_V for bes_df
    for filter_key in filter_keys:

        x = 1. / eff_wavelengths[filter_key]
        if x >= 1.1:
            y = x - 1.82
            a = 1 + y * (0.104 + y * (-0.609 + y * (0.701 + y * (
                        1.137 + y * (-1.718 + y * (-0.827 + y * (1.647 + y * (-0.505))))))))
            b = y * (1.952 + y * (2.908 + y * (-3.989 + y * (
                        -7.985 + y * (11.102 + y * (5.491 + y * (-10.805 + y * 3.347)))))))

        else:
            xpow = x ** 1.61
            a = 0.574 * xpow
            b = -0.527 * xpow

        Al_AV = a + b / R_V
        df[filter_key] += df['Av'] * Al_AV  # + 5*np.log10(bes_df['dist']*100)

    return df


def read_GSM_df(filename):
    df = pandas.read_csv('../outputfiles/' + filename)  # ,engine='python')
    return df


def read_BES_df(filename):
    names = ['Z087', 'Y106', 'J129', 'H158', 'W146', 'F184', 'mul', 'mub', 'Vr', 'UU', 'VV', 'WW',
        'Mv', 'CL', 'Typ', 'Teff', 'logg', 'Age', 'Mass', 'Mbol', 'Radius', '[Fe/H]', 'l(deg)',
        'b(deg)', 'RA2000.0', 'DEC2000.0', 'dist', 'x', 'y', 'z', 'Av', '[alpha/Fe]', 'weight']
    # names = ['J','Z087','H','W149','W169', 'mul', 'mub', 'Vr','UU', 'VV', 'WW','Mv' ,'CL', 'Typ', 'Teff', 'logg', 'Age','Mass','Mbol','Radius','[Fe/H]','l(deg)','b(deg)','RA2000.0','DEC2000.0','dist','x(kpc)','y','z' ,'Av','[alpha/Fe]' ]
    # df = pandas.read_csv('/home/stingray/johnson.7080/gulls/gulls_sj/lenses/EUCLID-lenses/'+filename,error_bad_lines=False,names=names,sep = '\s+')#,engine='python')
    df = pandas.read_csv('../outputfiles/' + filename, error_bad_lines=False, names=names,
        sep='\s+')  # ,engine='python')
    # Filter strings to NANs essentially
    for key in df.keys():
        df[key] = pandas.to_numeric(df[key], errors='coerce')
    # Drop rows with nans
    df = df.dropna()
    # reset index
    df.reset_index(drop=True)
    return df


def prep_plot_label(label):
    bad_chars = ['_', '(', '[', ']']

    for i in bad_chars:
        label = label.replace(i, '')

    return label


def plot_dist2(
        dict_list, xscale='linear', yscale='log', bins=None, dndbin=0, norm=1, plotfile=None,
        weight=1.523085e-07 / 0.000304617
        ):  # ,filter_bulge=True):
    """
    our model data frame and corresponding key
    gsm_df, gsm_key
    besancon dataframe and key
    bes_df, bes_key
    bins are None
    """

    # default will be distance
    # if bins is None or gsm_key is None or bes_df is None:
    #    #key for distance
    #     gsm_key = 'Dist'
    #    bes_key = 'dist'

    #    bins = np.linspace(0,20,100)

    # pdb.set_trace()
    # if weight:
    #    gsm_df['weight'] = np.ones(gsm_df.shape[0])*weight
    #    bes_df['weight'] = bes_df['weight']*weight
    # else:
    #    gsm_df['weight'] = np.ones(gsm_df.shape[0])
    #    bes_df['weight'] = np.ones(bes_df.shape[0])

    bin_centers = (bins[1:] + bins[:-1]) / 2
    hist_list = []

    fig, ax = plt.subplots()

    for df_dict in dict_list:
        df = df_dict[0]
        key = df_dict[1]
        label = df_dict[2]

        if weight and 'weight' not in df.keys():
            df['weight'] = np.ones(df.shape[0]) / weight
        # elif weight and 'weight' in df.keys():
        #    df['weight'] = df['weight']*weight
        # else:
        #    df['weight'] = np.ones(df.shape[0])

        if xscale == 'linear':
            tmp_hist, bins = np.histogram(df[key], bins=bins, weights=df['weight'])

        else:
            tmp_hist, bins = np.histogram(np.log10(df[key]), bins=bins, weights=df['weight'])

        ax.plot(bin_centers, tmp_hist, label=prep_plot_label(label))
        hist_list.append(tmp_hist)

    ax.set_yscale(yscale)
    ax.set_ylabel(r'N/deg$^2$')
    ax.set_xlabel(r'Histogram of Legend Values')
    plt.legend()
    if plotfile is not None:
        plt.savefig('./figures/' + plotfile)
    # plt.show()


def plot_dist(
        gsm_df, bes_df, gsm_key, bes_key, xscale='linear', yscale='log', bins=None, dndbin=0, norm=1
        ):
    """
    our model data frame and corresponding key
    gsm_df, gsm_key
    besancon dataframe and key
    bes_df, bes_key
    bins are None
    """

    # default will be distance
    if bins is None or gsm_key is None or bes_df is None:
        # key for distance
        gsm_key = 'Dist'
        bes_key = 'dist'

        bins = np.linspace(0, 20, 100)

    bin_centers = (bins[1:] + bins[:-1]) / 2
    if xscale == 'linear':
        gsm_hist, bins = np.histogram(gsm_df[gsm_key], bins=bins)
        bes_hist, bins = np.histogram(bes_df[bes_key], bins=bins)
    else:
        gsm_hist, bins = np.histogram(np.log10(gsm_df[gsm_key]), bins=bins)
        bes_hist, bins = np.histogram(np.log10(bes_df[bes_key]), bins=bins)

    fig, ax = plt.subplots()
    ax.plot(bin_centers, gsm_hist, label=prep_plot_label('GSM, ' + gsm_key))
    ax.plot(bin_centers, bes_hist, label=prep_plot_label('Bes, ' + bes_key))
    ax.set_yscale(yscale)
    ax.set_ylabel('N')
    ax.set_xlabel('Bins')
    plt.legend()
    # plt.show()


def plotmass(gsm_df, bes_df):
    plotfile = 'mass.pdf'
    dict_list = [[gsm_df, 'iMass', 'GSM, Initial Mass'],
        [gsm_df, 'Mass', 'GSM, Final Mass'],
        [bes_df, 'Mass', 'Bes, Final Mass']
        ]

    plot_dist2(dict_list, xscale='log', bins=np.linspace(-1.5, 1, 40), dndbin=0, norm=1,
        plotfile=plotfile)


def plotH(gsm_df, bes_df):
    plotfile = 'Hband.pdf'
    dict_list = [[gsm_df, '2MASS_H', 'GSM, H-band'],
        [bes_df, 'H', 'Bes, H-band']
        ]

    plot_dist2(dict_list, xscale='linear', bins=np.linspace(13, 40, 25), dndbin=0, norm=1,
        plotfile=plotfile)


def plotJ(gsm_df, bes_df):
    plotfile = 'Jband.pdf'
    dict_list = [[gsm_df, '2MASS_J', 'GSM, J-band'],
        [bes_df, 'J', 'Bes, J-band']
        ]

    plot_dist2(dict_list, xscale='linear', bins=np.linspace(13, 40, 25), dndbin=0, norm=1,
        plotfile=plotfile)


def plotCMD(gsm_df, bes_df):
    gsm_Z = gsm_df['Z087']
    gsm_J = gsm_df['2MASS_J']
    gsm_H = gsm_df['2MASS_H']
    # bes_J = bes_df['J']
    bes_Z = bes_df['Z087']
    bes_H = bes_df['H158']
    # pdb.set_trace()

    fig, ax = plt.subplots()
    # ax.plot(gsm_J-gsm_H,gsm_H,'.')
    # ax.plot(bes_J-bes_H,bes_H,'.')

    ax.plot(gsm_H - gsm_Z, gsm_H, '.')
    ax.plot(bes_H - bes_Z, bes_H, '.')

    ax.set_ylim(35, 8)
    ax.set_xlim(-8, 8)
    # plt.gca().invert_xaxis()

    ax.set_ylabel('H')
    ax.set_xlabel('H-Z')

    plt.savefig('./figures/cmd.pdf')


def plotAV(gsm_df, bes_df):
    # pdb.set_trace()
    AKs_AV = 0.13
    fig, ax = plt.subplots()
    av013 = np.array(gsm_df['A_Ks']) / 0.13
    av012 = np.array(gsm_df['A_Ks']) / 0.12

    ax.plot(gsm_df['Dist'], av013, '.', label='GSM, A(V), A(Ks)/A(V)=0.13')
    ax.plot(gsm_df['Dist'], av012, '.', label='GSM, A(V), A(Ks)/A(V)=0.12')
    # ax.plot(gsm_df['Dist'],gsm_df['A_Ks']/.12,'o',label='GSM, A(V), A(Ks)/A(V)=0.12')
    ax.plot(bes_df['dist'], bes_df['Av'], 'o', label='BES, A(V)')
    plt.legend()
    ax.set_ylabel('A(V)')
    ax.set_xlabel('Distance [kpc]')
    plt.savefig('plotAv.pdf')
    # plt.show()


def plotUU(gsm_df, bes_df):
    plotfile = 'UU.pdf'
    dict_list = [[gsm_df, 'UU', 'GSM, UU'],
        [bes_df, 'UU', 'Bes, UU']
        ]
    plot_dist2(dict_list, xscale='linear', bins=np.linspace(-100, 100, 60), dndbin=0, norm=1,
        plotfile=plotfile)


def plotVV(gsm_df, bes_df):
    plotfile = 'VV.pdf'
    dict_list = [[gsm_df, 'VV', 'GSM, VV'],
        [bes_df, 'VV', 'Bes, VV']
        ]
    plot_dist2(dict_list, xscale='linear', bins=np.linspace(-100, 100, 100), dndbin=0, norm=1,
        plotfile=plotfile)


def plotWW(gsm_df, bes_df):
    plotfile = 'WW.pdf'
    dict_list = [[gsm_df, 'WW', 'GSM, WW'],
        [bes_df, 'WW', 'Bes, WW']
        ]
    plot_dist2(dict_list, xscale='linear', bins=np.linspace(-100, 100, 100), dndbin=0, norm=1,
        plotfile=plotfile)


def plotdist(gsm_df, bes_df):
    plotfile = 'dist.pdf'
    dict_list = [[gsm_df, 'Dist', 'GSM, dist'],
        [bes_df, 'dist', 'Bes, dist']
        ]
    plot_dist2(dict_list, xscale='linear', bins=np.linspace(-1, 20, 100), dndbin=0, norm=1,
        plotfile=plotfile)


def plotX(gsm_df, bes_df):
    plotfile = 'X.pdf'
    dict_list = [[gsm_df, 'x', 'GSM, X'],
        [bes_df, 'x(kpc)', 'Bes, X']
        ]
    plot_dist2(dict_list, xscale='linear', bins=np.linspace(-15, 15, 60), dndbin=0, norm=1,
        plotfile=plotfile)


def plotY(gsm_df, bes_df):
    # pdb.set_trace()
    plotfile = 'Y.pdf'
    dict_list = [[gsm_df, 'y', 'GSM, Y'],
        [bes_df, 'y', 'Bes, Y']
        ]
    plot_dist2(dict_list, xscale='linear', bins=np.linspace(-0.1, 0, 60), dndbin=0, norm=1,
        plotfile=plotfile)


def plotZ(gsm_df, bes_df):
    plotfile = 'Z.pdf'
    dict_list = [[gsm_df, 'z', 'GSM, z'],
        [bes_df, 'z', 'Bes, z']
        ]
    plot_dist2(dict_list, xscale='linear', bins=np.linspace(-1, 1, 60), dndbin=0, norm=1,
        plotfile=plotfile)


def plotZ087(gsm_df, bes_df):
    plotfile = 'Z087.pdf'
    dict_list = [[gsm_df, 'Z087', 'GSM, Z-band'],
        [bes_df, 'Z087', 'Bes, Z-band']
        ]

    plot_dist2(dict_list, xscale='linear', bins=np.linspace(13, 40, 25), dndbin=0, norm=1,
        plotfile=plotfile)


def plotW146(gsm_df, bes_df):
    plotfile = 'W146.pdf'
    dict_list = [[gsm_df, 'W146', 'GSM, W146'],
        [bes_df, 'W146', 'Bes, W146']
        ]
    plot_dist2(dict_list, xscale='linear', bins=np.linspace(13, 40, 60), dndbin=0, norm=1,
        plotfile=plotfile)


def plot_bes_mags(bes_df):
    filter_keys = ['Z087', 'Y106', 'J129', 'H158', 'W146', 'F184']
    plotfile = 'bes_mags.pdf'
    dict_list = [[bes_df, 'Z087', 'Bes, Z087'],
        [bes_df, 'J', 'Bes, J-band'],
        [bes_df, 'Z087', 'Bes, Z-band'],
        [bes_df, 'H', 'Bes, H-band'],
        [bes_df, 'W149', 'Bes, W149'],
        [bes_df, 'W169', 'Bes, W169']
        ]
    plot_dist2(dict_list, xscale='linear', bins=np.linspace(13, 40, 25), dndbin=0, norm=1,
        plotfile=plotfile)


def plot_gsm_mags(bes_df):
    # kept as bes, should not matter
    plotfile = 'gsm_mags.pdf'
    dict_list = [[bes_df, '2MASS_J', 'Bes, J-band'],
        [bes_df, 'Z087', 'Bes, Z-band'],
        [bes_df, '2MASS_H', 'Bes, H-band'],
        [bes_df, 'W146', 'Bes, W146'],
        [bes_df, '2MASS_Ks', 'Bes, Ks-band']
        ]
    plot_dist2(dict_list, xscale='linear', bins=np.linspace(13, 40, 25), dndbin=0, norm=1,
        plotfile=plotfile)


def gsm_filter_pop(gsm_df, pop_id=10):
    gsm_df = gsm_df[gsm_df['pop'] < pop_id]
    return gsm_df


def bes_filter_pop(bes_df, pop_id=10):
    bes_df = bes_df[bes_df['Age'] < pop_id]
    return bes_df


def bes_filter_mass(bes_df, mass_lim=0.15):
    bes_df = bes_df[bes_df['Mass'] > mass_lim]

    return bes_df


def filter_dist(df, r=5, key='dist'):
    df = df[df[key] < r]
    return df


def dplot(gsm_df, bes_df):
    from mpl_toolkits import mplot3d
    gsm_skip = 30
    bes_skip = 30
    ax = plt.axes(projection='3d')
    ax.scatter3D(gsm_df['x'][::gsm_skip], gsm_df['y'][::gsm_skip], gsm_df['z'][::gsm_skip])
    ax.scatter3D(bes_df['x'][::bes_skip], bes_df['y'][::bes_skip], bes_df['z'][::bes_skip])
    plt.show()


if __name__ == '__main__':
    # pdb.set_trace()
    filename = '/lum_function_no_bulge/l-0.15b-3.2.csv'
    # filename = '/lum_function/l-0.15b-3.2.csv'
    gsm_df = read_GSM_df(filename)

    filename = '/lum_function_no_bulge/besfull-l-0.15b-3.2.wfirst'
    # filename = 'out-0014.wfirst'
    bes_df = read_BES_df(filename)
    pdb.set_trace()

    gsm_df = gsm_filter_pop(gsm_df, pop_id=8)
    bes_df = bes_filter_pop(bes_df, pop_id=8)
    dplot(gsm_df, bes_df)

    # for i in np.arange(10)+1:
    #    gdf = gsm_filter_pop(gsm_df,pop_id=i)
    #    bdf = bes_filter_pop(bes_df,pop_id=i)
    #    plotdist(gdf,bdf)
    #    plt.show()
    # bes_df = bes_filter_mass(bes_df)
    # gsm_df = filter_dist(gsm_df,r=4,key='Dist')
    # bes_df = filter_dist(bes_df,r=4)

    # plot_bes_mags(bes_df)
    # bes_df = apply_extinction(bes_df)
    plotdist(gsm_df, bes_df)
    plotmass(gsm_df, bes_df)

    plotW146(gsm_df, bes_df)
    plotZ087(gsm_df, bes_df)

    # plot_bes_mags(bes_df)
    # plot_gsm_mags(gsm_df)

    pdb.set_trace()

    plotdist(gsm_df, bes_df)
    plotAV(gsm_df, bes_df)
    plotCMD(gsm_df, bes_df)
    plotdist(gsm_df, bes_df)
    # plotAV(gsm_df,bes_df)
    plotH(gsm_df, bes_df)
    plotX(gsm_df, bes_df)
    plotY(gsm_df, bes_df)
    plotZ(gsm_df, bes_df)
    plotJ(gsm_df, bes_df)
    plotH(gsm_df, bes_df)
    # plotJ(gsm_df,bes_df)
    plotZband(gsm_df, bes_df)
    plotUU(gsm_df, bes_df)
    plotVV(gsm_df, bes_df)
    plotWW(gsm_df, bes_df)

    pdb.set_trace()
