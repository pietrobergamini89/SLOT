import numpy as np
import holoviews as hv
import param
from holoviews import opts
import pandas as pd
from bokeh.models import Div
from holoviews import streams

from holoviews.operation.datashader import rasterize
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
from astropy.wcs import WCS
from astropy.table import Table

import panel as pn
import io

import os
from subprocess import Popen, PIPE
import string

from bokeh.models import HoverTool

import glob
from astropy import units as u

import zipfile

import pickle

import time

pn.extension(loading_spinner='dots', loading_color='#00aa41')
hv.extension('bokeh')
pn.param.ParamMethod.loading_indicator = True


def clusters(dir_in='./CLUSTERS'):
    cluster_dict = dict()

    for cl in os.listdir(dir_in):

        model_dict = dict()

        if '.' not in cl:
            cluster_dict[cl] = dict()
            image_dict = dict()

            for images_path in os.listdir(dir_in + '/' + cl + '/IMAGES'):
                if '.fits' in images_path:
                    image_dict[images_path] = str(dir_in + '/' + cl + '/IMAGES/'+images_path)

            cluster_dict[cl]['IMAGES'] = image_dict

            for models_path in os.listdir(dir_in + '/' + cl + '/LENS_MODELS'):

                if '.' not in models_path:

                    model_dict[models_path] = dict()
                    with open(dir_in + '/' + cl + '/LENS_MODELS/' + models_path + '/SLOT.txt') as conf:
                        lines = conf.readlines()
                        for line in lines:
                            if ':' in line:
                                if 'REF_COORDS' in line:
                                    model_dict[models_path][line.split(': ')[0]] = \
                                        np.asarray([line.split(': ')[1].replace('\n', '').split(',')[0],
                                                    line.split(': ')[1].replace('\n', '').split(',')[1]], dtype=float)
                                else:
                                    model_dict[models_path][line.split(': ')[0]] = line.split(': ')[1].replace('\n', '')

                cluster_dict[cl]['LENS_MODELS'] = model_dict

    return cluster_dict


if not os.path.isdir('./computations'):
    os.mkdir('./computations')


# def image_preparator(dir_in='./CLUSTERS/IMAGES'):
#     im_dict = dict()
#     images = os.listdir(dir_in)
#     print(images)
#     for i, image in enumerate(images):
#         path = dir_in+'/'+image
#         try:
#             with fits.open(path) as hdu:
#                 im_dict[image] = [np.asarray(hdu[0].data, dtype=float), hdu[0].header]
#         except:
#             pass
#     return im_dict


def coord_diff(ra_ref, dec_ref, ra, dec):
    x = -(ra - ra_ref) * np.cos(dec_ref / 180.0 * np.pi) * 3600.0
    y = (dec - dec_ref) * 3600.0
    return x, y


# def lensing_preparator(dir_in='./CLUSTERS/MACS0416/LENS_MODELS'):
#     mod_dict = dict()
#
#     for models_path in os.listdir(dir_in):
#         if '.' not in models_path:
#             mod_dict[models_path] = dict()
#             with open(dir_in + '/' + models_path + '/SLOT.txt') as conf:
#                 lines = conf.readlines()
#                 for line in lines:
#                     if ':' in line:
#                         if 'REF_COORDS' in line:
#                             mod_dict[models_path][line.split(': ')[0]] = \
#                                 np.asarray([line.split(': ')[1].replace('\n', '').split(',')[0],
#                                             line.split(': ')[1].replace('\n', '').split(',')[1]], dtype=float)
#                         else:
#                             mod_dict[models_path][line.split(': ')[0]] = line.split(': ')[1].replace('\n', '')
#
#     return mod_dict

# start_time = time.time()
# print("--- %s seconds ---" % (time.time() - start_time))

if os.path.isfile('bk_files.pkl'):
    clusters_dic = pickle.load(open('bk_files.pkl', 'rb'))
else:
    clusters_dic = clusters(dir_in='./CLUSTERS')

    with open('bk_files.pkl', 'wb') as f:
        pickle.dump(clusters_dic, f)



# current_cluster = ['']
current_cluster = [list(clusters_dic.keys())[0]]
wcs = ['']
# current_band = ['']
current_band = [list(clusters_dic[current_cluster[0]]['IMAGES'].keys())[0]]
# current_model = ['']
current_model = [list(clusters_dic[current_cluster[0]]['LENS_MODELS'].keys())[0]]

# class CLUSTERselector(param.Parameterized):
#
#     cluster = param.ObjectSelector(default=list(clusters_dic.keys())[0], objects=list(clusters_dic.keys()), label='')
#
#     @param.depends('cluster')
#     def change_cluster(self):
#         current_cluster[0] = self.cluster


# class MODELselector(param.Parameterized):
#
#     model = param.ObjectSelector(default=list(clusters_dic[current_cluster[0]]['LENS_MODELS'].keys())[0],
#                                  objects=list(clusters_dic[current_cluster[0]]['LENS_MODELS'].keys()), label='')
#
#     @param.depends('model')
#     def change_model(self):
#         current_model[0] = self.model


class MAINplot(param.Parameterized):

    cluster = param.ObjectSelector(default=list(clusters_dic.keys())[0], objects=list(clusters_dic.keys()), label='')

    model = param.ObjectSelector(default=list(clusters_dic[current_cluster[0]]['LENS_MODELS'].keys())[0],
                                 objects=list(clusters_dic[current_cluster[0]]['LENS_MODELS'].keys()), label='')

    filter_selector = param.ObjectSelector(default=list(clusters_dic[current_cluster[0]]['IMAGES'].keys())[0],
                                           objects=list(clusters_dic[current_cluster[0]]['IMAGES'].keys()),
                                           label='Filter')

    scale_selector = param.Range(default=(0, 0.1), bounds=(0, 1), step=0.01, label='Scale')

    references = param.ListSelector(objects=["Create region"])

    df = pd.DataFrame({'ID': np.asarray([], dtype=str), 'z': [], 'X': [], 'Y': [], 'RA': [], 'DEC': []})

    df_results = pd.DataFrame({'ID': np.asarray([], dtype=str), 'z': [], 'X_best': [], 'Y_best': [], 'RA_best': [],
                               'dRA_16': [], 'dRA_50': [], 'dRA_84': [], 'DEC_best': [], 'dDEC_16': [], 'dDEC_50': [],
                               'dDEC_84': [], 'Xsource_best': [], 'Ysource_best': [], 'RAsource_best': [],
                               'DECsource_best': [], 'AMP_best': [], 'AMP_16': [], 'AMP_50': [], 'AMP_84': [],
                               'DTIME_best': [], 'DTIME_16': [], 'DTIME_50': [], 'DTIME_84': [],
                               'PARITY': np.asarray([], dtype=str), 'SEP_16': [], 'SEP_50': [], 'SEP_84': []})

    im = param.Action(lambda x: x.param.trigger('im'))
    gal = param.ListSelector(objects=['Cluster galaxies'])

    clear_df = param.Action(lambda x: x.param.trigger('clear_df'))

    id_input = param.String(doc='ID')
    z_input = param.String(doc='z')
    ra_input = param.String(doc='RA')
    dec_input = param.String(doc='Dec')

    insert = param.Action(lambda x: x.param.trigger('insert'))

    update_ref = param.Action(lambda x: x.param.trigger('update_ref'))

    mcmc = param.Boolean(False, doc="MCMC")

    compute_values = param.Action(lambda x: x.param.trigger('compute_values'))

    input_table = param.Parameter()

    table_file_name = param.String(default="Input_table.csv", doc='The filename to save to.')

    table_results_file_name = param.String(default="Table_results.csv", doc='The filename to save to.')

    maps_selector = param.ObjectSelector(default='Magnification', objects=['Magnification', 'Convergence', 'Shear',
                                                                   'Surface mass density', 'Deflection components'],
                                         label='Map type:')

    z_map = param.String(doc='z')
    pix_map = param.String(doc='Pixels')
    ra_map = param.String(doc='RA')
    dec_map = param.String(doc='Dec')

    size_map = param.Number(100, bounds=(1, 200), step=0.5)

    map_name = param.String(default="Input_table.csv", doc='The filename to save to.')

    generate_map = param.Action(lambda x: x.param.trigger('generate_map'))

    xc = 0
    yc = 0

    def __init__(self, **params):
        super().__init__(**params)
        self.save_table = pn.widgets.FileDownload(name="", filename=self.table_file_name, callback=self.save_table_file,
                                                  margin=(0, 0, 25, 5), button_type="primary")

        self.save_rtable = pn.widgets.FileDownload(name="", filename=self.table_results_file_name,
                                                   callback=self.save_rtable_file, margin=(0, 0, 15, 5),
                                                   button_type="primary")

        with fits.open(clusters_dic[current_cluster[0]]['IMAGES'][current_band[0]]) as hdu:
            self.data = np.asarray(hdu[0].data, dtype=float)
            wcs[0] = WCS(hdu[0].header)

    def inizializer(self, len1, len2):
        a = np.empty((len1, len2))
        a[:] = np.nan
        return a

    def notunique(self, a):
        unique = []
        repeted = []
        for ele in a:
            if ele not in unique:
                unique.append(ele)
            else:
                repeted.append(ele)
        return repeted

    def bayesImageResults(self, path, dist_best):

        obs = Table.read(dist_best, format='ascii', header_start=1)
        ra_ref = float(obs.meta['comments'][0].split(' ')[2])
        dec_ref = float(obs.meta['comments'][0].split(' ')[3])
        obs['X'] = ra_ref - obs['X'] / (np.cos(dec_ref / 180.0 * np.pi) * 3600.0)
        obs['Y'] = obs['Y'] / 3600.0 + dec_ref

        dist = glob.glob(path + '/dist*')
        dist.sort()
        cc = np.zeros(len(obs['ID']))

        res_AMP = self.inizializer(len(dist), len(cc))
        res_DTIME = self.inizializer(len(dist), len(cc))
        res_sep = self.inizializer(len(dist), len(cc))
        res_x = self.inizializer(len(dist), len(cc))
        res_y = self.inizializer(len(dist), len(cc))

        num = np.asarray(obs['ID'])

        for i, f in enumerate(dist):

            t = Table.read(f, format='ascii', header_start=1)

            ra_match = ra_ref - t['X'] / (np.cos(dec_ref / 180.0 * np.pi) * 3600.0)
            dec_match = t['Y'] / 3600.0 + dec_ref

            for n, number in enumerate(np.unique(num)):
                for par in ['++', '+-', '-+', '--']:

                    mask_dist = (t['ID'] == number) * (t['PARITY'] == par)
                    mask_obs = (num == number) * (obs['PARITY'] == par)

                    if len(obs[mask_obs]) != 0 and len(t[mask_dist]) != 0:

                        c = SkyCoord(ra=obs['X'][mask_obs] * u.degree, dec=obs['Y'][mask_obs] * u.degree)

                        catalog = SkyCoord(ra=ra_match[mask_dist] * u.degree, dec=dec_match[mask_dist] * u.degree)

                        idx, sep, sep3 = c.match_to_catalog_sky(catalog)
                        sep = sep.arcsec

                        if not np.array_equal(num[mask_obs],
                                              np.asarray(t['ID'][mask_dist][idx], dtype=float)) or not \
                                ~np.array_equal(obs['PARITY'][mask_obs], np.asarray(t['PARITY'][mask_dist][idx],
                                                                                    dtype=str)):
                            quit()

                        res_AMP[i][mask_obs] = np.asarray(t['AMP'][mask_dist][idx], dtype=float)
                        res_DTIME[i][mask_obs] = np.asarray(t['DTIME'][mask_dist][idx], dtype=float)
                        res_sep[i][mask_obs] = np.asarray(sep, dtype=float)
                        res_x[i][mask_obs] = np.asarray(t['X'][mask_dist][idx], dtype=float)
                        res_y[i][mask_obs] = np.asarray(t['Y'][mask_dist][idx], dtype=float)

                        if len(np.unique(idx)) != len(num[mask_obs]):
                            mask_Fnunique = np.zeros(len(idx), dtype=bool)

                            for uni in self.notunique(idx):
                                mask_notunique = (idx == uni)
                                min_sep = (sep == np.min(sep[mask_notunique]))
                                mask_Fnunique[mask_notunique * (~min_sep)] = True

                            mask_obs[mask_obs] = mask_Fnunique

                            res_AMP[i][mask_obs] = np.nan
                            res_DTIME[i][mask_obs] = np.nan
                            res_sep[i][mask_obs] = np.nan
                            res_x[i][mask_obs] = np.nan
                            res_y[i][mask_obs] = np.nan

        percentiles_AMP = np.zeros((len(res_AMP[0]), 3))
        percentiles_DTIME = np.zeros((len(res_AMP[0]), 3))
        percentiles_sep = np.zeros((len(res_AMP[0]), 3))
        x = np.zeros((len(res_AMP[0]), 3))
        y = np.zeros((len(res_AMP[0]), 3))

        for im, amp in enumerate(res_AMP[0]):

            mask_nan = ~np.ma.masked_invalid(res_x[:, im]).mask
            try:
                percentiles_AMP[im] = np.percentile(res_AMP[:, im][mask_nan], [16, 50, 84])
                percentiles_DTIME[im] = np.percentile(res_DTIME[:, im][mask_nan], [16, 50, 84])
                percentiles_sep[im] = np.percentile(res_sep[:, im][mask_nan], [16, 50, 84])
                x[im] = np.percentile(res_x[:, im][mask_nan], [16, 50, 84])
                y[im] = np.percentile(res_y[:, im][mask_nan], [16, 50, 84])

            except:
                percentiles_AMP[im] = np.nan
                percentiles_DTIME[im] = np.nan
                percentiles_sep[im] = np.nan
                x[im] = np.nan
                y[im] = np.nan

        return percentiles_AMP, percentiles_DTIME, percentiles_sep, x, y

    def doubletap(self, x, y):
        if x is not None and y is not None:
            self.xc = x
            self.yc = y
            radec = pixel_to_skycoord(x, y, wcs[0])
            self.ra_input = str(radec.ra.deg)
            self.dec_input = str(radec.dec.deg)
            self.ra_map = str(radec.ra.deg)
            self.dec_map = str(radec.dec.deg)

        return hv.Points([(self.xc, self.yc)]).options({'Points': {'marker': '+', 'size': 10, 'color': 'red'}})

    @param.depends('cluster', watch=True)
    def change_cluster(self):
        current_cluster[0] = self.cluster
        self.df = pd.DataFrame({'ID': np.asarray([], dtype=str), 'z': [], 'X': [], 'Y': [], 'RA': [], 'DEC': []})

        self.df_results = pd.DataFrame({'ID': np.asarray([], dtype=str), 'z': [], 'X_best': [], 'Y_best': [],
                                        'RA_best': [], 'dRA_16': [], 'dRA_50': [], 'dRA_84': [], 'DEC_best': [],
                                        'dDEC_16': [], 'dDEC_50': [], 'dDEC_84': [], 'Xsource_best': [],
                                        'Ysource_best': [], 'RAsource_best': [], 'DECsource_best': [], 'AMP_best': [],
                                        'AMP_16': [], 'AMP_50': [], 'AMP_84': [], 'DTIME_best': [], 'DTIME_16': [],
                                        'DTIME_50': [], 'DTIME_84': [], 'PARITY': np.asarray([], dtype=str),
                                        'SEP_16': [], 'SEP_50': [], 'SEP_84': []})
        self.param.model.objects = list(clusters_dic[self.cluster]['LENS_MODELS'].keys())
        self.model = list(clusters_dic[self.cluster]['LENS_MODELS'].keys())[0]
        self.param.filter_selector.objects = list(clusters_dic[self.cluster]['IMAGES'].keys())
        self.filter_selector = list(clusters_dic[self.cluster]['IMAGES'].keys())[0]

        current_band[0] = list(clusters_dic[self.cluster]['IMAGES'].keys())[0]
        current_model[0] = list(clusters_dic[self.cluster]['LENS_MODELS'].keys())[0]



    @param.depends('model')
    def change_model(self):
        current_model[0] = self.model

    # @param.depends('filter_selector', 'scale_selector', 'references', 'df', 'df_results', 'ra_map', 'dec_map')
    @param.depends('filter_selector', 'scale_selector', 'references', 'update_ref', 'gal', 'im', 'insert',
                   'compute_values', 'input_table', 'clear_df')
    def image_creator(self):
        try:
            if self.filter_selector != current_band[0]:
                with fits.open(clusters_dic[self.cluster]['IMAGES'][self.filter_selector]) as hdu:
                    self.data = np.asarray(hdu[0].data, dtype=float)
                    wcs[0] = WCS(hdu[0].header)
                    current_band[0] = self.filter_selector

            x_ax = np.arange(wcs[0].pixel_shape[0])
            y_ax = np.arange(wcs[0].pixel_shape[1])

            # MAIN IMAGE

            mimage = hv.Image((x_ax, y_ax, self.data))
            image = rasterize(mimage.opts(aspect='equal', responsive=True, cmap='Greys', clim=(self.scale_selector[0],
                                                                                               self.scale_selector[1]),
                                          colorbar=True, tools=['zoom_in', 'zoom_out', 'undo', 'redo', 'reset']),
                              precompute=True)

            double_tap = hv.streams.DoubleTap(source=mimage)
            dtap = hv.DynamicMap(self.doubletap, streams=[double_tap])

            # INPUT MULTIPLE IMAGES

            hover_im = HoverTool(tooltips=[('z', '@z')])
            input_images = hv.Points(self.df, kdims=['X', 'Y']).opts(color='limegreen', marker='+', size=10, tools=[hover_im])
            input_images_labels = hv.Labels((self.df['X'], self.df['Y'], self.df['ID'])).opts(opts.Labels(
                text_color='limegreen', text_font_size='12pt'))

            # OUTPUT MULTIPLE IMAGES

            hover_res = HoverTool(tooltips=[('AMP', '@AMP_best'), ('DTIME', '@DTIME_best'), ('PARITY', '@PARITY')])

            output_images = hv.Points(self.df_results, kdims=['X_best', 'Y_best']).opts(color='orange', marker='+', size=10,
                                                                                        tools=[hover_res])
            output_images_labels = hv.Labels((self.df_results['X_best'], self.df_results['Y_best'],
                                              self.df_results['ID'])).opts(opts.Labels(text_color='orange',
                                                                                       text_font_size='12pt'))

            sources = hv.Points(self.df_results, kdims=['Xsource_best', 'Ysource_best']).opts(color='blue', marker='+',
                                                                                              size=10)
            sources_labels = hv.Labels((self.df_results['Xsource_best'], self.df_results['Ysource_best'],
                                        self.df_results['ID'])).opts(opts.Labels(text_color='blue', text_font_size='12pt'))

            # INPUT CLUSTER GALAXIES

            if self.gal is not None:
                if 'Cluster galaxies' in self.gal:
                    df_gal = Table.read('CLUSTERS/'+current_cluster[0]+'/LENS_MODELS/'+current_model[0]+'/' +
                                        clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['CLUSTER_GALAXIES'],
                                        format='ascii', names=['ID', 'RA', 'DEC', 'a', 'b', 't', 'm_F160', 'l']).to_pandas()

                    coords_cg = SkyCoord(df_gal['RA'], df_gal['DEC'], unit='deg', frame='fk5')
                    xypix_cg = skycoord_to_pixel(coords_cg, wcs[0])

                    df_gal['X'] = xypix_cg[0]
                    df_gal['Y'] = xypix_cg[1]
                    hover_gal = HoverTool(tooltips=[('ID', '@ID'), ('Mag_F160W', '@m_F160')])
                    output_gal = hv.Scatter(df_gal, kdims=['X', 'Y']).opts(line_color='magenta', color='none', marker='o',
                                                                           size=20, tools=[hover_gal])
                else:
                    output_gal = hv.Points([[0, 0]]).opts(color='white', marker='+', size=0)
            else:
                output_gal = hv.Points([[0, 0]]).opts(color='white', marker='+', size=0)

            if self.references is not None:
                if 'Create region' in self.references:
                    if self.ra_map != '' and self.dec_map != '':
                        pixscale = wcs[0].pixel_scale_matrix[1, 1] * 3600

                        cref = skycoord_to_pixel(SkyCoord(clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['REF_COORDS'][0],
                                                          clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['REF_COORDS'][1], unit='deg',
                                                          frame='fk5'), wcs[0])

                        x_arcsec, y_arcsec = coord_diff(clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['REF_COORDS'][0],
                                                        clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['REF_COORDS'][1], float(self.ra_map),
                                                        float(self.dec_map))

                        xmapmin = cref[0] + x_arcsec / pixscale - (self.size_map / 2) / pixscale
                        xmapmax = cref[0] + x_arcsec / pixscale + (self.size_map / 2) / pixscale
                        ymapmin = cref[1] + y_arcsec / pixscale - (self.size_map / 2) / pixscale
                        ymapmax = cref[1] + y_arcsec / pixscale + (self.size_map / 2) / pixscale

                        map = hv.Curve([(xmapmin, ymapmin), (xmapmin, ymapmax), (xmapmax, ymapmax), (xmapmax, ymapmin),
                                        (xmapmin, ymapmin)]).opts(color='lightblue')

                        return image * map * input_images * input_images_labels * dtap * output_gal * output_images * \
                            output_images_labels * sources * sources_labels

            return image * input_images * input_images_labels * dtap * output_gal * output_images * output_images_labels * \
                sources * sources_labels
        except:
            return image

    @param.depends('filter_selector', 'scale_selector')
    def histo_scale(self):
        if self.filter_selector != current_band[0]:
            with fits.open(clusters_dic[self.cluster]['IMAGES'][self.filter_selector]) as hdu:
                self.data = np.asarray(hdu[0].data, dtype=float)
                wcs[0] = WCS(hdu[0].header)
            # self.data = image_dict[self.filter_selector][0]
            current_band[0] = self.filter_selector

        histo, edges = np.histogram(self.data, range=[self.scale_selector[0], self.scale_selector[1]], density=False,
                                    bins=500)

        return hv.Histogram((edges, np.log10(histo))).opts(axiswise=True, height=170, width=297, toolbar='above',
                                                           xlabel='pix', ylabel='10Ey')

    @param.depends('im', watch=True)
    def images_model(self):
        table_imlt_astropy = Table.read('CLUSTERS/'+current_cluster[0]+'/LENS_MODELS/'+current_model[0]+'/' +
                                clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['MULTIPLE_IMAGES'],
                                format='ascii', names=['ID', 'RA', 'DEC', 'a', 'b', 't', 'z', 'l'])
        table_imlt_astropy.keep_columns(['ID', 'RA', 'DEC', 'z'])
        table_imlt = table_imlt_astropy.to_pandas()

        coords_imlt = SkyCoord(table_imlt['RA'], table_imlt['DEC'], unit='deg', frame='fk5')
        xypix_imlt = skycoord_to_pixel(coords_imlt, wcs[0])
        table_imlt['X'] = xypix_imlt[0]
        table_imlt['Y'] = xypix_imlt[1]
        table_imlt = table_imlt[['ID', 'z', 'X', 'Y', 'RA', 'DEC']]

        self.df = self.df.append(table_imlt)

    @param.depends('df')
    def view_df(self):
        return pn.widgets.DataFrame(self.df, width=1000, height=170, show_index=False, auto_edit=True, name='',
                                    widths={'ID': 40, 'z': 60, 'X': 120, 'Y': 120, 'RA': 120, 'DEC': 120})

    @param.depends('df_results')
    def view_df_results(self):
        return pn.widgets.DataFrame(self.df_results, width=1550, height=170, autosize_mode='none', show_index=False,
                                    auto_edit=False, name='', widths={'ID': 40, 'z': 60, 'X_best': 120, 'Y_best': 120,
                                                                      'RA_best': 120, 'dRA_16': 120, 'dRA_50': 120,
                                                                      'dRA_84': 120, 'DEC_best': 120, 'dDEC_16': 120,
                                                                      'dDEC_50': 120, 'dDEC_84': 120,
                                                                      'Xsource_best': 120, 'Ysource_best': 120,
                                                                      'RAsource_best': 120, 'DECsource_best': 120,
                                                                      'AMP_best': 120, 'AMP_16': 120, 'AMP_50': 120,
                                                                      'AMP_84': 120, 'DTIME_best': 120, 'DTIME_16': 120,
                                                                      'DTIME_50': 120, 'DTIME_84': 120, 'PARITY': 60,
                                                                      'SEP_16': 120, 'SEP_50': 120, 'SEP_84': 120})

    @param.depends('clear_df', watch=True)
    def clear(self):
        self.df = pd.DataFrame({'ID': np.asarray([], dtype=str), 'z': [], 'X': [], 'Y': [], 'RA': [], 'DEC': []})

        self.df_results = pd.DataFrame({'ID': np.asarray([], dtype=str), 'z': [], 'X_best': [], 'Y_best': [],
                                        'RA_best': [], 'dRA_16': [], 'dRA_50': [], 'dRA_84': [], 'DEC_best': [],
                                        'dDEC_16': [], 'dDEC_50': [], 'dDEC_84': [], 'Xsource_best': [],
                                        'Ysource_best': [], 'RAsource_best': [], 'DECsource_best': [], 'AMP_best': [],
                                        'AMP_16': [], 'AMP_50': [], 'AMP_84': [], 'DTIME_best': [], 'DTIME_16': [],
                                        'DTIME_50': [], 'DTIME_84': [], 'PARITY': np.asarray([], dtype=str),
                                        'SEP_16': [], 'SEP_50': [], 'SEP_84': []})

    @param.depends('insert', watch=True)
    def add_value(self):
        try:
            if self.xc != 0 or self.yc != 0:
                if ',' in self.z_input:
                    val = self.z_input.split(',')
                    z_multiple = np.arange(float(val[0]), float(val[1]), float(val[2]))
                    id_multiple_list = []

                    for i in np.arange(len(z_multiple)):
                        id_multiple_list.append(self.id_input[:-1] + str(i+1) + self.id_input[-1])

                    add_df = pd.DataFrame({'ID': np.asarray(id_multiple_list, dtype=str),
                                           'z': z_multiple,
                                           'X': np.full(len(z_multiple), float(self.xc)),
                                           'Y': np.full(len(z_multiple), float(self.yc)),
                                           'RA': np.full(len(z_multiple), float(self.ra_input)),
                                           'DEC': np.full(len(z_multiple), float(self.dec_input))})

                    self.df = pd.concat([self.df, add_df], ignore_index=True)

                else:
                    self.df = self.df.append({'ID': self.id_input,
                                              'z': float(self.z_input),
                                              'X': float(self.xc),
                                              'Y': float(self.yc),
                                              'RA': float(self.ra_input),
                                              'DEC': float(self.dec_input)}, ignore_index=True)
        except:
            pass

    @param.depends('input_table', watch=True)
    def read_file(self):
        try:
            table_im_astropy = Table.read(io.BytesIO(self.input_table), format='ascii')
            table_im_astropy.keep_columns(['ID', 'RA', 'DEC', 'z'])
            table_im = table_im_astropy.to_pandas()
            coords_im = SkyCoord(table_im['RA'], table_im['DEC'], unit='deg', frame='fk5')
            xypix_im = skycoord_to_pixel(coords_im, wcs[0])
            table_im['X'] = xypix_im[0]
            table_im['Y'] = xypix_im[1]
            table_im = table_im[['ID', 'z', 'X', 'Y', 'RA', 'DEC']]
            self.df = self.df.append(table_im)
        except:
            pass

    @param.depends('compute_values', watch=True)
    def compute(self):
        try:
            self.df_results = pd.DataFrame({'ID': np.asarray([], dtype=str), 'z': [], 'X_best': [], 'Y_best': [],
                                            'RA_best': [], 'dRA_16': [], 'dRA_50': [], 'dRA_84': [], 'DEC_best': [],
                                            'dDEC_16': [], 'dDEC_50': [], 'dDEC_84': [], 'Xsource_best': [],
                                            'Ysource_best': [], 'RAsource_best': [], 'DECsource_best': [], 'AMP_best': [],
                                            'AMP_16': [], 'AMP_50': [], 'AMP_84': [], 'DTIME_best': [], 'DTIME_16': [],
                                            'DTIME_50': [], 'DTIME_84': [], 'PARITY': np.asarray([], dtype=str),
                                            'SEP_16': [], 'SEP_50': [], 'SEP_84': []})

            im = Table()
            im['col1'] = np.array(self.df['ID'], dtype=str)
            im['col2'] = np.array(self.df['RA'], dtype=float)
            im['col3'] = np.array(self.df['DEC'], dtype=float)
            im['col4'] = np.full(len(self.df['RA']), 0.1, dtype=float)
            im['col5'] = np.full(len(self.df['RA']), 0.1, dtype=float)
            im['col6'] = np.full(len(self.df['RA']), 0, dtype=float)
            im['col7'] = np.array(self.df['z'], dtype=float)
            im['col8'] = np.full(len(self.df['RA']), 25, dtype=float)

            dirname = 'computations/' + str(np.random.randint(1, 1000000)) + '/'
            os.system('mkdir ' + dirname[:-1])

            os.system('cp ' + 'CLUSTERS/'+current_cluster[0]+'/LENS_MODELS/'+current_model[0]+'/' + clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['BAYES_BEST']
                      + ' ' + dirname + 'bayes.dat')

            os.system('cp ' + 'CLUSTERS/'+current_cluster[0]+'/LENS_MODELS/'+current_model[0]+'/' +clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['INPUT']
                      + ' ' + dirname + clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['INPUT'])

            os.system('cp ' + 'CLUSTERS/'+current_cluster[0]+'/LENS_MODELS/'+current_model[0]+'/' +clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['CLUSTER_GALAXIES']
                      + ' ' + dirname + clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['CLUSTER_GALAXIES'])

            im.meta['comments'] = ['REFERENCE 0']
            im.write(dirname + clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['MULTIPLE_IMAGES'], format='ascii.fixed_width_no_header',
                     delimiter='\t', overwrite=True)

            bayesimage = Popen(['bayesImage ' + clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['INPUT']], stdout=PIPE, stderr=PIPE,
                               shell=True, cwd=dirname)
            bayesimage.wait()

            sources = Table.read(dirname + 'images/source0.dat', format='ascii')
            dist = Table.read(dirname + 'images/dist0.dat', format='ascii', header_start=1)
            ra_ref = float(dist.meta['comments'][0].split(' ')[2])
            dec_ref = float(dist.meta['comments'][0].split(' ')[3])

            x_b = np.asarray(dist['X'])
            y_b = np.asarray(dist['Y'])

            dist['X'] = ra_ref - dist['X'] / (np.cos(dec_ref / 180.0 * np.pi) * 3600.0)
            dist['Y'] = dist['Y'] / 3600.0 + dec_ref

            sources['col2'] = ra_ref - sources['col2'] / (np.cos(dec_ref / 180.0 * np.pi) * 3600.0)
            sources['col3'] = sources['col3'] / 3600.0 + dec_ref

            names = dist['ID'].astype(str)
            redshift = np.zeros(len(names))
            ras = np.zeros(len(names))
            decs = np.zeros(len(names))

            l = list(string.ascii_lowercase)
            idp = -1
            il = 0

            for i, ids in enumerate(dist['ID']):
                mask_sources = (sources['col1'] == ids)
                redshift[i] = float(sources['col7'][mask_sources])
                ras[i] = float(sources['col2'][mask_sources])
                decs[i] = float(sources['col3'][mask_sources])

                if ids != idp:
                    idp = ids
                    il = 0

                names[i] = str(ids) + l[il]

                il += 1

            dist.replace_column('ID', names)

            coords = SkyCoord(dist['X'], dist['Y'], unit='deg', frame='fk5')
            xypix = skycoord_to_pixel(coords, wcs[0])

            coordss = SkyCoord(ras, decs, unit='deg', frame='fk5')
            xypixs = skycoord_to_pixel(coordss, wcs[0])

            if self.mcmc:

                dirname2 = 'computations/' + str(np.random.randint(1, 1000000)) + '/'

                os.system('mkdir ' + dirname2[:-1])

                os.system('cp '+'CLUSTERS/'+current_cluster[0]+'/LENS_MODELS/'+current_model[0]+'/' + clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['BAYES_RANDOM']
                          + ' ' + dirname2 + 'bayes.dat')

                os.system('cp '+'CLUSTERS/'+current_cluster[0]+'/LENS_MODELS/'+current_model[0]+'/' + clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['INPUT']
                          + ' ' + dirname2 + clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['INPUT'])

                os.system('cp '+'CLUSTERS/'+current_cluster[0]+'/LENS_MODELS/'+current_model[0]+'/' + clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['CLUSTER_GALAXIES']
                          + ' ' + dirname2 + clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['CLUSTER_GALAXIES'])

                im.write(dirname2 + clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['MULTIPLE_IMAGES'], format='ascii.fixed_width_no_header',
                         delimiter='\t', overwrite=True)

                bayesimage = Popen(['bayesImage ' + clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['INPUT']], stdout=PIPE, stderr=PIPE,
                                   shell=True, cwd=dirname2)
                bayesimage.wait()

                percentiles_AMP, percentiles_DTIME, percentiles_sep, x, y = self.bayesImageResults(dirname2 + '/images',
                                                                                              dirname + '/images/dist0.dat')

                df3 = pd.DataFrame(
                    {'ID': dist['ID'], 'z': redshift, 'X_best': xypix[0], 'Y_best': xypix[1], 'RA_best': dist['X'],
                     'dRA_16': np.round(x[:, 0], 4) - x_b, 'dRA_50': np.round(x[:, 1], 4) - x_b,
                     'dRA_84': np.round(x[:, 2], 4) - x_b, 'DEC_best': dist['Y'], 'dDEC_16': np.round(y[:, 0], 4) - y_b,
                     'dDEC_50': np.round(y[:, 1], 4) - y_b, 'dDEC_84': np.round(y[:, 2], 4) - y_b,
                     'Xsource_best': xypixs[0], 'Ysource_best': xypixs[1], 'RAsource_best': ras, 'DECsource_best': decs,
                     'AMP_best': dist['AMP'], 'AMP_16': np.round(percentiles_AMP[:, 0], 3),
                     'AMP_50': np.round(percentiles_AMP[:, 1], 3), 'AMP_84': np.round(percentiles_AMP[:, 2], 3),
                     'DTIME_best': dist['DTIME'], 'DTIME_16': np.round(percentiles_DTIME[:, 0], 3),
                     'DTIME_50': np.round(percentiles_DTIME[:, 1], 3), 'DTIME_84': np.round(percentiles_DTIME[:, 2], 3),
                     'PARITY': dist['PARITY'], 'SEP_16': np.round(percentiles_sep[:, 0], 3),
                     'SEP_50': np.round(percentiles_sep[:, 1], 3), 'SEP_84': np.round(percentiles_sep[:, 2], 3)})

                os.system('rm -r ' + dirname2[:-1])

            else:
                df3 = pd.DataFrame(
                    {'ID': dist['ID'], 'z': redshift, 'X_best': xypix[0], 'Y_best': xypix[1], 'RA_best': dist['X'],
                     'DEC_best': dist['Y'], 'Xsource_best': xypixs[0], 'Ysource_best': xypixs[1], 'RAsource_best': ras,
                     'DECsource_best': decs, 'AMP_best': dist['AMP'], 'DTIME_best': dist['DTIME'],
                     'PARITY': dist['PARITY']})

            os.system('rm -r ' + dirname[:-1])

            self.df_results = df3
        except:
            pass

    @param.depends('table_file_name', watch=True)
    def update_filename_table(self):
        self.save_table.filename = self.table_file_name

    @param.depends('save_table')
    def save_table_file(self):
        self.save_table.filename = self.table_file_name
        sio = io.StringIO()
        self.df.to_csv(sio, index=False)
        sio.seek(0)
        return sio

    @param.depends('table_results_file_name', watch=True)
    def update_filename_rtable(self):
        self.save_rtable.filename = self.table_results_file_name

    @param.depends('save_rtable')
    def save_rtable_file(self):
        self.save_rtable.filename = self.table_results_file_name
        sio = io.StringIO()
        self.df_results.to_csv(sio, index=False)
        sio.seek(0)
        return sio

    def table_download_widget(self):
        return self.save_table

    def rtable_download_widget(self):
        return self.save_rtable

    def map_widget(self):
        return self.mapd

    @param.depends('maps_selector', watch=True)
    def map_sel(self):
        if self.maps_selector == 'Surface mass density':
            self.z_map = clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['Z']

    @param.depends('generate_map')
    def map_generator(self):
        try:
            if self.z_map != '' and self.pix_map != '' and self.ra_map != '' and self.dec_map != '':

                command = ''
                if self.maps_selector == 'Magnification':
                    command += '    ampli 1 ' + self.pix_map + ' ' + self.z_map + ' Magnification_z' + self.z_map + '_Npix' + \
                               self.pix_map + '.fits \n'
                    name = 'Magnification_z' + self.z_map + '_Npix' + self.pix_map + '.fits'
                # if self.maps_selector == 'Magnification modulus':
                #     command += '    ampli 2 ' + self.pix_map + ' ' + self.z_map + ' AbsAmpli_z' + self.z_map + '_Npix' + \
                #                self.pix_map + '.fits \n'
                #     name = 'AbsAmpli_z' + self.z_map + '_Npix' + self.pix_map + '.fits'
                if self.maps_selector == 'Convergence':
                    command += '    ampli 5 ' + self.pix_map + ' ' + self.z_map + ' Convergence_Npix' + self.pix_map + \
                               '.fits \n'
                    name = 'Convergence_Npix' + self.pix_map + '.fits'
                if self.maps_selector == 'Shear':
                    command += '    ampli 6 ' + self.pix_map + ' ' + self.z_map + ' Shear_Npix' + self.pix_map + '.fits \n'
                    name = 'Shear_Npix' + self.pix_map + '.fits'
                if self.maps_selector == 'Surface mass density':
                    self.z_map = clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['Z']
                    command += '    mass 4 ' + self.pix_map + ' ' + str(float(clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['Z'])) + ' Mass_z' + str(float(clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['Z'])) + '_Npix' + \
                               self.pix_map + '.fits \n'
                    name = 'Mass_z' + str(float(clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['Z'])) + '_Npix' + self.pix_map + '.fits'

                # if self.maps_selector == 'Mass_pix2':
                #     command += '    mass 3 ' + self.pix_map + ' ' + '0.396' + ' Mass_pix2_z' + '0.396' + '_Npix' + \
                #                self.pix_map + '.fits \n'
                #     name = 'Mass_pix2_z' + self.z_map + '_Npix' + self.pix_map + '.fits'
                if self.maps_selector == 'Deflection components':
                    command += '    dpl 1 ' + self.pix_map + ' ' + self.z_map + ' dx_z' + self.z_map + '_Npix' + \
                               self.pix_map + '.fits' + ' dy_z' + self.z_map + '_Npix' + self.pix_map + '.fits \n'
                    name = 'dx_z' + self.z_map + '_Npix' + self.pix_map + '.fits' + ',dy_z' + self.z_map + '_Npix' + \
                           self.pix_map + '.fits'

                dirname_maps = 'computations/' + str(np.random.randint(1, 1000000)) + '/'

                os.system('mkdir ' + dirname_maps[:-1])

                with open('CLUSTERS/'+current_cluster[0]+'/LENS_MODELS/'+current_model[0]+'/' + clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['BEST'], 'r') as f:
                    lines = f.readlines()
                lines.insert(14, command)

                x_arcsec, y_arcsec = coord_diff(clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['REF_COORDS'][0],
                                                clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['REF_COORDS'][1], float(self.ra_map),
                                                float(self.dec_map))

                xmapmin = np.round(x_arcsec - float(self.size_map) / 2, 1)
                xmapmax = np.round(x_arcsec + float(self.size_map) / 2, 1)
                ymapmin = np.round(y_arcsec - float(self.size_map) / 2, 1)
                ymapmax = np.round(y_arcsec + float(self.size_map) / 2, 1)

                lines[-6] = '    xmin     ' + str(xmapmin) + '\n'
                lines[-5] = '    xmax     ' + str(xmapmax) + '\n'
                lines[-4] = '    ymin     ' + str(ymapmin) + '\n'
                lines[-3] = '    ymax     ' + str(ymapmax) + '\n'

                with open(dirname_maps + 'best.par', 'w') as g:
                    lines = "".join(lines)
                    g.write(lines)

                ltmap = Popen(['lenstool best.par -n'], stdout=PIPE, stderr=PIPE, shell=True, cwd=dirname_maps)
                ltmap.wait()

                if ',' in name:
                    with zipfile.ZipFile(dirname_maps + 'deflections.zip', 'w') as zipObj:
                        zipObj.write(dirname_maps + name.split(',')[0])
                        zipObj.write(dirname_maps + name.split(',')[1])
                    file_name = dirname_maps + 'deflections.zip'
                    name = 'deflections.zip'
                else:
                    file_name = dirname_maps + name

                self.mapd = pn.widgets.FileDownload(file_name, embed=True, filename=name, button_type="primary", width=285,
                                                    margin=(5, 0, 10, 11))

                os.system('rm -r ' + dirname_maps[:-1])

                return self.mapd
        except:
            pass
        try:
            os.system('rm -r ' + dirname_maps[:-1])
        except:
            pass


class LensFiles(param.Parameterized):

    def __init__(self, **params):
        super().__init__(**params)

        self.file_download_par = pn.widgets.FileDownload('CLUSTERS/'+current_cluster[0]+'/LENS_MODELS/'+current_model[0]+'/' + clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['INPUT'],
                                                         embed=True,
                                                         filename=clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['INPUT'],
                                                         button_type="primary", width=285)

        self.file_download_images = pn.widgets.FileDownload('CLUSTERS/'+current_cluster[0]+'/LENS_MODELS/'+current_model[0]+'/' + clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['MULTIPLE_IMAGES'],
                                                            embed=True,
                                                            filename=clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['MULTIPLE_IMAGES'],
                                                            button_type="primary",
                                                            width=285)

        self.file_download_members = pn.widgets.FileDownload('CLUSTERS/'+current_cluster[0]+'/LENS_MODELS/'+current_model[0]+'/' + clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['CLUSTER_GALAXIES'],
                                                             embed=True,
                                                             filename=clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['CLUSTER_GALAXIES'],
                                                             button_type="primary",
                                                             width=285)

        self.file_download_bayes = pn.widgets.FileDownload('CLUSTERS/'+current_cluster[0]+'/LENS_MODELS/'+current_model[0]+'/' + clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['BAYES'],
                                                           embed=True,
                                                           filename=clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['BAYES'],
                                                           button_type="primary",
                                                           width=285)

        self.file_download_chires = pn.widgets.FileDownload('CLUSTERS/'+current_cluster[0]+'/LENS_MODELS/'+current_model[0]+'/' + clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['CHIRES'],
                                                            embed=True,
                                                            filename=clusters_dic[current_cluster[0]]['LENS_MODELS'][current_model[0]]['CHIRES'],
                                                            button_type="primary",
                                                            width=285)

    def par_widget(self):
        return self.file_download_par

    def images_widget(self):
        return self.file_download_images

    def members_widget(self):
        return self.file_download_members

    def bayes_widget(self):
        return self.file_download_bayes

    def chires_widget(self):
        return self.file_download_chires



mainplot = MAINplot()
lens_files = LensFiles()
# model_selector = MODELselector()

title = Div(text="<b>Cluster:</b>",
                 style={'font-size': '16pt', 'color': 'black'}
                 )

subtitle = Div(text="<b>Lens Model:</b>",
               style={'font-size': '12pt', 'color': 'black'}
               )

input_files = Div(text="<b>Input files</b>",
                  style={'font-size': '14pt', 'color': 'black'}
                  )

output_files = Div(text="<b>Best-fit lens model</b>",
                   style={'font-size': '14pt', 'color': 'black'}
                   )

credits = Div(text="<b>Info:</b> pietro.bergamini@unimi.it",
               style={'font-size': '12pt', 'color': 'black'}
               )

fixed_widgets = pn.WidgetBox(
                        pn.Param(mainplot, parameters=["filter_selector"], show_name=False),
                        mainplot.histo_scale,
                        pn.Param(mainplot, widgets={"scale_selector": {"widget_type": pn.widgets.EditableRangeSlider, 'width': 280}}, parameters=["scale_selector"], show_name=False),
                        pn.layout.Divider(),
                        pn.Param(mainplot, parameters=["gal"],
                                 widgets={"gal": {"widget_type": pn.widgets.CheckButtonGroup}}, show_name=False),
                        pn.layout.Divider(),
                            )


tab1 = pn.WidgetBox(
                 pn.layout.Divider(),
                 pn.Param(mainplot, parameters=["im"], widgets={"im": {"widget_type": pn.widgets.Button,
                                                                       'name': 'Model images'}}, show_name=False),
                 pn.layout.Divider(),
                 pn.Row(
                        pn.Param(mainplot, widgets={"id_input": {'widget_type': pn.widgets.TextInput,
                                                                 'placeholder': 'ID', 'name': ''}},
                                 parameters=["id_input"], width=145, show_name=False),
                        pn.Param(mainplot, widgets={"z_input": {'widget_type': pn.widgets.TextInput,
                                                                'placeholder': 'z', 'name': ''}},
                                 parameters=["z_input"], width=145, show_name=False)
                       ),
                 pn.Row(
                        pn.Param(mainplot, widgets={"ra_input": {'widget_type': pn.widgets.TextInput,
                                                                 'placeholder': 'RA', 'name': ''}},
                                 parameters=["ra_input"], width=145, show_name=False),
                        pn.Param(mainplot, widgets={"dec_input": {'widget_type': pn.widgets.TextInput,
                                                                  'placeholder': 'Dec', 'name': ''}},
                                 parameters=["dec_input"], width=145, show_name=False)
                       ),
                 pn.Param(mainplot.param, widgets={"insert": {'widget_type': pn.widgets.Button,
                                                              "button_type": "success", 'loading_indicator': True}},
                          parameters=["insert"], show_name=False),
                 pn.Param(mainplot.param, widgets={"clear_df": {'widget_type': pn.widgets.Button,
                                                                "button_type": "danger", 'name': 'Clear'}},
                          parameters=["clear_df"], show_name=False),
                 pn.layout.Divider(),
                 pn.Param(mainplot, widgets={"input_table": pn.widgets.FileInput}, parameters=["input_table"],
                          show_name=False),
                 pn.layout.Divider(),
                 pn.Param(mainplot, widgets={"mcmc": {'widget_type': pn.widgets.Checkbox, 'name': 'Compute errors'}},
                          parameters=["mcmc"], show_name=False),
                 pn.Param(mainplot, widgets={"compute_values": {'widget_type': pn.widgets.Button,
                                                                "button_type": "warning", 'name': 'Compute',
                                                                'margin': (5, 0, 10, 6), 'width': 285}},
                          parameters=["compute_values"], show_name=False)
                    )

tab2 = pn.WidgetBox(pn.layout.Divider(), input_files, lens_files.par_widget, lens_files.members_widget,
                    lens_files.images_widget, pn.layout.Divider(), output_files, lens_files.bayes_widget,
                    lens_files.chires_widget, pn.layout.Divider())

tab3 = pn.WidgetBox(
                   pn.layout.Divider(),
                   pn.Param(mainplot, parameters=["maps_selector"], show_name=False),
                   pn.Param(mainplot, widgets={"references": {'widget_type': pn.widgets.CheckButtonGroup,
                                                              'width': 280, 'margin': (10, 0, 1, 9)}},
                            parameters=["references"], show_name=False),
                   pn.Param(mainplot, widgets={"update_ref": {'widget_type': pn.widgets.Button,
                                                              "button_type": "success", 'name': 'Update region',
                                                              'margin': (0, 0, 10, 9), 'width': 280}},
                            parameters=["update_ref"], show_name=False),
                   pn.Row(
                         pn.Param(mainplot, widgets={"z_map": {'widget_type': pn.widgets.TextInput,
                                                               'placeholder': 'z', 'name': ''}}, parameters=["z_map"],
                                  show_name=False, width=145),
                         pn.Param(mainplot, widgets={"pix_map": {'widget_type': pn.widgets.TextInput,
                                                                 'placeholder': 'Npix', 'name': ''}},
                                  parameters=["pix_map"], show_name=False, width=145)
                         ),
                   pn.Row(
                         pn.Param(mainplot, widgets={"ra_map": {'widget_type': pn.widgets.TextInput,
                                                                'placeholder': 'RA', 'name': ''}},
                                  parameters=["ra_map"], width=145, show_name=False),
                         pn.Param(mainplot, widgets={"dec_map": {'widget_type': pn.widgets.TextInput,
                                                                 'placeholder': 'Dec', 'name': ''}},
                                  parameters=["dec_map"], width=145, show_name=False)
                         ),
                   pn.Param(mainplot, widgets={"size_map": {"widget_type": pn.widgets.EditableFloatSlider, 'name': 'Map side (arcsec)',
                                                            'width': 280}}, parameters=["size_map"], show_name=False),
                   pn.Param(mainplot, widgets={"generate_map": {'widget_type': pn.widgets.Button,
                                                                "button_type": "warning", 'width': 285,
                                                                'margin': (15, 0, 10, 6)}},
                            parameters=["generate_map"], show_name=False),
                   mainplot.map_generator
                   )

widgets = pn.Column(title,
                    pn.Param(mainplot,
                             widgets={"cluster": {'widget_type': pn.widgets.Select}},
                             parameters=["cluster"], show_name=False),
                    subtitle,
                    pn.Param(mainplot,
                             widgets={"model": {'widget_type': pn.widgets.Select}},
                             parameters=["model"], show_name=False),
                    pn.layout.Divider(margin=(0, 30, 0, 5)),
                    fixed_widgets,
                    pn.Tabs(('Images', tab1), ('Model', tab2), ('Maps', tab3), width=350),
                    )

pn.Column(pn.Row(widgets,
                 mainplot.image_creator),
          mainplot.view_df,
          pn.layout.Divider(margin=(1, 0, 0, 0)),
          pn.Param(mainplot, widgets={"table_file_name": {"widget_type": pn.widgets.TextInput, 'name': '',
                                                          'margin': (1, 0, 0, 0)}}, parameters=["table_file_name"],
                   show_name=False, width=145),
          mainplot.table_download_widget,
          mainplot.view_df_results,
          pn.layout.Divider(margin=(1, 0, 0, 0)),
          pn.Param(mainplot, widgets={"table_results_file_name": {"widget_type": pn.widgets.TextInput, 'name': '',
                                                                  'margin': (1, 0, 0, 0)}},
                   parameters=["table_results_file_name"], show_name=False, width=145),
          mainplot.rtable_download_widget,
          height=200, width=200, sizing_mode='scale_width', loading_indicator=True
          ).servable()
