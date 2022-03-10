#
# Mobilion Acorn REST API services
# Copyright (c) 2022 MOBILion Systems, Inc.
# Author: Harry Collins
#
"""Test CCS calibration components & endpoint"""
import os
import json
from typing import List, Dict, Tuple, Optional
import matplotlib.pyplot as plt
from pdfme import build_pdf
from canonical_names import CanonicalName
from data_transfer_classes import ChannelValue
import Polarity
from Polarity import Polarity
from mbifile_wrapper import MBIFileManager
from mobiliondata.datareader import DataReader, FrameCollection
from models2 import QuartileMethod, PeakPicker, GaussianFitter
from analyte import Analyte
from channelvalueserializer import ChannelValueSerializer


class MBIfile():
    """test reading the MBI file and looking at it's data"""
    analyte_data = [{"name": "118", "mz_value": 118.086255, "polarity": "positive", "ccs": 121.3, "ccs_range": 0.2},
                    {"name": "322", "mz_value": 322.048121, "polarity": "positive", "ccs": 153.75, "ccs_range": 0.23},
                    {"name": "622", "mz_value": 622.028960, "polarity": "positive", "ccs": 202.96, "ccs_range": 0.27},
                    {"name": "922", "mz_value": 922.009798, "polarity": "positive", "ccs": 243.64, "ccs_range": 0.3},
                    {"name": "1222", "mz_value": 1221.990637, "polarity": "positive", "ccs": 282.2, "ccs_range": 0.47},
                    {"name": "1522", "mz_value": 1521.971475, "polarity": "positive", "ccs": 316.96, "ccs_range": 0.6},
                    {"name": "1822", "mz_value": 1821.952313, "polarity": "positive", "ccs": 351.25, "ccs_range": 0.62},
                    {"name": "2122", "mz_value": 2121.933152, "polarity": "positive", "ccs": 383.03, "ccs_range": 0.64},
                    {"name": "2422", "mz_value": 2421.913990, "polarity": "positive", "ccs": 412.96, "ccs_range": 0.58},
                    {"name": "2722", "mz_value": 2721.894829, "polarity": "positive", "ccs": 441.21, "ccs_range": 0.59},
                    {"name": "113", "mz_value": 112.985587, "polarity": "negative", "ccs": 108.23, "ccs_range": 0.2},
                    {"name": "302", "mz_value": 301.998139, "polarity": "negative", "ccs": 140.04, "ccs_range": 0.29},
                    {"name": "602", "mz_value": 601.978977, "polarity": "negative", "ccs": 180.77, "ccs_range": 0.21},
                    {"name": "1034", "mz_value": 1033.988109, "polarity": "negative", "ccs": 255.34, "ccs_range": 0.32},
                    {"name": "1334", "mz_value": 1333.968947, "polarity": "negative", "ccs": 284.76, "ccs_range": 0.31},
                    {"name": "1634", "mz_value": 1633.949786, "polarity": "negative", "ccs": 319.03, "ccs_range": 0.7},
                    {"name": "1934", "mz_value": 1933.930624, "polarity": "negative", "ccs": 352.55, "ccs_range": 0.27},
                    {"name": "2234", "mz_value": 2233.911463, "polarity": "negative", "ccs": 380.74, "ccs_range": 0.31},
                    {"name": "2534", "mz_value": 2533.892301, "polarity": "negative", "ccs": 412.99, "ccs_range": 0.31},
                    {"name": "2834", "mz_value": 2833.873139, "polarity": "negative", "ccs": 432.62, "ccs_range": 0.35}]

    @staticmethod
    def readbacks_to_channels(readbacks: str, use_db: bool = False) -> List[ChannelValue]:
        """given a list of readbacks (potentially legacy data) turn it into channel values"""
        readback_channels: dict = json.loads(readbacks)
        channels: List[ChannelValue] = []
        if isinstance(readback_channels, list) and 'channel' in readback_channels[0]:  # this is the old format
            if not use_db:  # using SimpleTestCase, so we have no DB access, just assume +
                sg_inputs = [item for item in readback_channels if item['label'] == 'SLIM Guard']
                if sg_inputs:
                    sg_readback = [item for item in sg_inputs[0]['inputs'] if item['label'] == 'readback']
                    if sg_readback:
                        slim_guard = float(sg_readback[0]['value'])
                        channels.append(ChannelValue(CanonicalName.QTOF_POLARITY.value,
                                                     value=Polarity.POSITIVE.value if slim_guard >= 0.0 else Polarity.NEGATIVE.value))
                sw_inputs = [item for item in readback_channels if item['label'] == 'Separation Traveling Wave']
                if sw_inputs:
                    sw_frequency = [item for item in sw_inputs[0]['inputs'] if item['label'] == 'frequency'][0]['value']
                    channels.append(ChannelValue(canonical_name='Separation Traveling Wave.frequency', value=str(sw_frequency)))
        else:
            if not use_db:
                if isinstance(readback_channels, list):
                    channels = [ChannelValue(item['canonical_name'], item['value']) for item in readback_channels]
                else:
                    channels = [ChannelValue(item, readback_channels[item]) for item in readback_channels]
            else:
                serializer = ChannelValueSerializer(data=readback_channels, many=True)
                serializer.is_valid(raise_exception=True)
                channels = serializer.validated_data

        channels_dict: Dict[str, ChannelValue] = {item.canonical_name: item for item in channels}
        if CanonicalName.QTOF_POLARITY.value not in channels_dict:
            # create a polarity
            polarity = Polarity.POSITIVE.value
            try:
                if 'SLIM Guard.readback' in channels_dict:
                    if float(channels_dict['SLIM Guard.readback'].value) < 0.0:
                        polarity = polarity.NEGATIVE.value
                elif 'SLIM Guard.TWBias' in channels_dict:
                    if float(channels_dict['SLIM Guard.TWBias'].value) < 0.0:
                        polarity = Polarity.NEGATIVE.value
            except KeyError:
                pass

            if channels:
                channels.append(ChannelValue(CanonicalName.QTOF_POLARITY, value=polarity, timestamp=channels[0].timestamp))

        return channels

    def read_frames(self, mbi_data):
        """test reading in frames for stock data"""
        # avoid the DB so this is fast, create our analytes (just for the polarity of the test)
        all_experiments_data = []
        all_analytes: List[Analyte] = []
        for item in MBIfile.analyte_data:
            analyte = Analyte(name=item['name'], mz_value=item['mz_value'], polarity=item['polarity'], ccs=item['ccs'], ccs_range=item['ccs_range'])
            all_analytes.append(analyte)


        mbi_files: {} = mbi_data
        for mbi_file_path in mbi_files:
            plot_name: str = [item for item in mbi_files[mbi_file_path]][0]
            expected_data: dict = mbi_files[mbi_file_path][plot_name]
            try:
                parts: List[str] = mbi_file_path.split('\\')
            except:
                file_path = mbi_file_path
                file_name = mbi_file_path
            else:
                file_path = mbi_file_path
                name = len(parts) - 1
                file_name: str = parts[name]
            header: str = f"\nCalibration file -> {file_name}"
            print(f"{header}\n{'='*len(header)}")
            with MBIFileManager(file_path, mode='r') as mbi_file:
                data_reader: DataReader = DataReader(mbi_file)
                global_dict = data_reader.global_metadata()
                ccs_cal = json.loads(global_dict.get('cal-ccs'))
                at_surfing = round(ccs_cal.get('at_surfing'), 3)
                ccs_min = ccs_cal.get('min')
                ccs_max = ccs_cal.get('max')
                ccs_coefficients = ccs_cal.get('coefficients')
                analytes = ccs_cal.get('peaks')

                frames: FrameCollection = data_reader.frame_collection(1)
                tof_pulse_period: float = float(frames[0].metadata('frm-dt-period'))  # length of a TOF pulse frame in seconds
                # the 'bin_width' is defined as the digitizer period, so in this case about 0.5ns

                readbacks: List[ChannelValue] = self.readbacks_to_channels(data_reader.global_metadata(key='readbacks'))

                # sum the frames, the 'drift_times' are the arrival times
                frames: FrameCollection = data_reader.frame_collection(1, -1)  # get all the frames

                # there are some values that we can pull from the method stashed in the MBI file
                method_channels = self.readbacks_to_channels(readbacks=frames[0].metadata('frm-method-state'))
                qtof_polarity: ChannelValue = [item for item in method_channels if item.canonical_name == CanonicalName.QTOF_POLARITY][0]
                sw_frequency: float = float([item for item in method_channels if item.canonical_name == 'Separation Traveling Wave.frequency'][0].value)
                at_swf: float = 1000.0 * (13.0 / (sw_frequency * 0.009))  # ms

                # now that we know the polarity, filter the analytes
                analytes = [analyte for analyte in all_analytes if analyte.polarity == qtof_polarity.value]

                drift_times: dict = {}
                for analyte in analytes:
                    dt_for_analyte = frames.dt_time(mz_start=analyte.mz_value-0.5, mz_stop=analyte.mz_value+0.5)
                    drift_times[analyte.name] = {item: dt_for_analyte[item] for item in dt_for_analyte if dt_for_analyte[item]}

                    if len(drift_times[analyte.name]) < 3 or analyte.name in ['322']:
                        print(f"{plot_name}: {analyte.name} missing data!")
                        continue

                    pp = PeakPicker(drift_times[analyte.name], sample_period=tof_pulse_period/10.0, min_peaks=3)
                    at_peak: float = pp.peak()
                    try:
                        print(f"{analyte.name}, expected/measured AT (ms) = {expected_data[analyte.name]['at']}/{round(at_peak, 3)}")
                    except KeyError:
                        pass

                    # let's plot the arrival times

                    try:
                        expected_at: float = expected_data[analyte.name]['at']
                    except KeyError:
                        pass
                    expected_gauss_centroid: float = expected_data[analyte.name]['gauss']
                    expected_centroid: float = expected_data[analyte.name]['centroid']

                    qm = QuartileMethod(drift_times[analyte.name], peak=at_peak)
                    pp = PeakPicker(qm.eligible_samples, sample_period=tof_pulse_period/10.0, min_peaks=3)
                    oversampled: {} = pp.oversampled(sample_period=tof_pulse_period/10.0)

                    x = [item for item in oversampled]
                    y = [oversampled[item] for item in x]
                    gauss_fitter: GaussianFitter = GaussianFitter(x=x, y=y, center=at_peak, amplitude=pp.max_intensity, sigma=qm.median - qm.Q1)
                    result = gauss_fitter.fit()
                    peak_center: float = gauss_fitter.center
                    peak_fwhm: float = gauss_fitter.fwhm
                    print(result.fit_report())

                    # using the newly calculate FWHM, let's reduce the list of eligible samples further for display purposes
                    if peak_center < 0.0:
                        print(f"Gaussian for {analyte.name} has a negative centroid!")
                        continue

                    x = [item for item in oversampled if item <= peak_center + peak_fwhm*2.0]
                    if not x:
                        print(f"No data for {file_name}, {analyte.sample}, center/fwhm = {peak_center}/{peak_fwhm}")
                        continue

                    y = [oversampled[item] for item in x]
                    # re-fit to the newly constrained data
                    gauss_fitter: GaussianFitter = GaussianFitter(x=x, y=y, center=at_peak, amplitude=pp.max_intensity, sigma=qm.median - qm.Q1)
                    result = gauss_fitter.fit()
                    print(result.fit_report())
                    save_gauss_report = result.fit_report()

                    plt.plot(x, y, c='blue')
                    plt.xlabel('arrival time (ms)')
                    plt.ylabel('intensity')
                    plt.title(f"{analyte.name}, {file_name} {'(+)' if qtof_polarity.value == 'positive' else '(-)'}")

                    y_min, y_max = plt.gca().get_ylim()
                    x_min, x_max = plt.gca().get_xlim()

                    legend_plot = []
                    legend_text = []
                    peak_legend,  = plt.plot([at_peak, at_peak], [y_min, y_max], c='red', linestyle='dashed')
                    legend_plot.append(peak_legend)
                    legend_text.append(f'peak {round(at_peak, 3)}')
                    if expected_at:
                        exp_legend, = plt.plot([expected_at, expected_at], [y_min, y_max], c='green', linestyle='dashed')
                        legend_plot.append(exp_legend)
                        legend_text.append(f'expected {round(expected_at, 3)}')
                    if expected_gauss_centroid:
                        gauss_legend, = plt.plot([expected_gauss_centroid, expected_gauss_centroid], [y_min, y_max], c='black', linestyle='dashdot')
                        legend_plot.append(gauss_legend)
                        legend_text.append(f'gauss {round(expected_gauss_centroid, 3)} vs {round(peak_center,3)}')
                    else:
                        gauss_legend, = plt.plot([peak_center, peak_center], [y_min, y_max], c='black', linestyle='dashdot')
                        legend_plot.append(gauss_legend)
                        legend_text.append(f'gauss {round(peak_center,3)} (calculated)')

                    if expected_centroid:
                        centroid_legend, = plt.plot([expected_centroid, expected_centroid], [y_min, y_max], c='pink', linestyle='dotted')
                        legend_plot.append(centroid_legend)
                        legend_text.append(f'centroid {round(expected_centroid, 3)}')

                    plt.legend(legend_plot, legend_text)  # display legends we have

                    fwhm_lower, fwhm_upper = pp.FWHM(peak_at=at_peak)
                    plt.plot([fwhm_lower[0], fwhm_lower[0]], [0.0, fwhm_lower[1]], c='black', linestyle='solid')
                    plt.plot([fwhm_upper[0], fwhm_upper[0]], [0.0, fwhm_upper[1]], c='black', linestyle='solid')
                    plt.plot(x, result.best_fit, c='orange')

                    if peak_center <= at_swf * 1.1:
                        plt.text((x_min + x_max) / 2.0, (y_min + y_max) / 2.0, 'surfing',
                                 fontsize=40, color='gray', alpha=0.5,
                                 ha='center', va='center', rotation='30')
                    plt.savefig(fr'graphs\{analyte.name}_{file_name}.jpg')
                    plt.close()
                    graph = {
                        "style": {
                            "margin_bottom": 5, "text_align": "j",
                            "page_size": "letter", "margin": [60, 50]
                        },
                        "formats": {
                            "url": {"c": "blue", "u": 1},
                            "title": {"b": 1, "s": 13}
                        },
                        "running_sections": {
                            "header": {
                                "x": "left", "y": 20, "height": "top", "style": {"text_align": "r"},
                                "content": [{".b": "This is a header"}]
                            },
                            "footer": {
                                "x": "left", "y": 740, "height": "bottom", "style": {"text_align": "c"},

                            }
                        },
                        "sections": [
                            {
                                "style": {"page_numbering_style": "roman"},
                                "running_sections": ["footer"],
                                "content": [

                                    {
                                        ".": " ", "style": "title", "label": "title1",
                                        "outline": {"level": 1, "text": "none"}
                                    },

                                    {"image": fr'graphs\{analyte.name}_{file_name}.jpg'},

                                ]

                            },

                        ]
                    }
                    report = {
                        "style": {
                            "margin_bottom": 5, "text_align": "j",
                            "page_size": "letter", "margin": [60, 50]
                        },
                        "formats": {
                            "url": {"c": "blue", "u": 1},
                            "title": {"b": 1, "s": 13}
                        },
                        "running_sections": {
                            "header": {
                                "x": "left", "y": 20, "height": "top", "style": {"text_align": "r"},
                                "content": [{".b": "This is a header"}]
                            },
                            "footer": {
                                "x": "left", "y": 740, "height": "bottom", "style": {"text_align": "c"},

                            }
                        },
                        "sections": [
                            {
                                "style": {"page_numbering_style": "roman"},
                                "running_sections": ["footer"],
                                "content": [

                                    {
                                        ".": "CE07_Andys_CCS_Test_that_WILL_work.mbi", "style": "title", "label": "title1",
                                        "outline": {"level": 1, "text": "wutthisdo"}
                                    },

                                    f'Report for {analyte.name}_{file_name}\n \nCCS: {analyte.ccs}\n{save_gauss_report}'
                                ]

                            },

                        ]
                    }

                    with open(fr'pdfs\{analyte.name}_{file_name}_graph.pdf', 'wb') as f:
                        build_pdf(graph, f)
                    with open(fr'pdfs\{analyte.name}_{file_name}_report.pdf', 'wb') as f:
                        build_pdf(report, f)
                    temp_coordinate_dict = {}
                    for idx in range(len(x)):
                        temp_coordinate_dict.update({x[idx]: y[idx]})
            document = {
                "style": {
                    "margin_bottom": 5, "text_align": "j",
                    "page_size": "letter", "margin": [60, 50]
                },
                "formats": {
                    "url": {"c": "blue", "u": 1},
                    "title": {"b": 1, "s": 13}
                },
                "running_sections": {
                    "header": {
                        "x": "left", "y": 20, "height": "top", "style": {"text_align": "r"},
                        "content": [{".b": "This is a header"}]
                    },
                    "footer": {
                        "x": "left", "y": 740, "height": "bottom", "style": {"text_align": "c"},

                    }
                },
                "sections": [
                    {
                        "style": {"page_numbering_style": "roman"},
                        "running_sections": ["footer"],
                        "content": [

                            {
                                ".": "CE07_Andys_CCS_Test_that_WILL_work.mbi", "style": "title", "label": "title1",
                                "outline": {"level": 1, "text": "none"}
                            },

                            f'at_surfing: {at_surfing}ms\nCCS-Min: {ccs_min}\nCCS-Max: {ccs_max}\nCoefficients: {ccs_coefficients}'
                        ]

                    },

                ]
            }
            with open('endpage.pdf',
                      'wb') as f:
                build_pdf(document, f)
