from typing import Dict, Any
import json
import matplotlib.pyplot as plt
from mobiliondata import mbi
from mobiliondata.defs import DataCubeKeys

target_list = ["Mass Flow.gas type", "Mass Flow.setpoint", "Mass Flow.pressure", "Funnel Pressure.pressure-serial", "SLIM Bias.DCBias", "SLIM Bias.readback", "Separation Traveling Wave.frequency", "Separation Traveling Wave.amplitude", "Separation Traveling Wave.readback", "On Board Accumulation Traveling Wave.amplitude", "On Board Accumulation Traveling Wave.readback", "Short Path Gate.control", "Short Path Gate.gate", "Short Path Gate.readback", "Long Path Gate.control", "Long Path Gate.gate", "Long Path Gate.readback", "Funnel In.DCBias", "Funnel In.readback", "Funnel Exit.DCBias", "Funnel Exit.readback", "Funnel CL.DCBias", "Funnel CL.readback", "Deflector.TWBias", "Deflector.readback", "Exit CL.DCBias", "Exit CL.readback", "Quad Bias.DCBias", "Quad Bias.readback", "SLIM Top RF.drive", "SLIM Top RF.rf+", "SLIM Top RF.rf-", "SLIM Bottom RF.drive", "SLIM Bottom RF.rf+", "SLIM Bottom RF.rf-", "SLIM Guard.TWBias", "SLIM Guard.readback", "Quad RF.drive", "Quad RF.rf+", "Quad RF.frequency", "Funnel RF.drive", "Funnel RF.rf+", "Funnel RF.rf-", "Funnel RF.frequency", "QTOF Polarity.DIO", "Timing Table.release", "Timing Table.trap", "Timing Table.acquire", "Timing Table.accumulations", "SLIM.path_length"]
mbifile = mbi.Open('mbidata\CE07_Andys_CCS_Test_that_WILL_work.mbi')
data_reader = mbi.DataReader(mbifile)
frame_col = data_reader.frame_collection(1, -1)
session = frame_col[15]
frame_range = (0, 47)
ADC_value = list(session.index_counts[:])
gate_start_stop = list(session.index_positions[:])
unpacked_dict = {}
frames = [frame_col[i] for i in range(len(frame_col))]
max_scans = 0
max_idx = 0
scan_entries_total = []
num_gates_total = []
num_samples_total = []
unpacked_total = []
scan_tic_total = []
adc_scans_total = []
gate_scans_total = []
for frm_idx, frm in enumerate(frames):
    if len(frm) > max_scans:
        max_scans = len(frm)
        max_idx = frm_idx
    session = frame_col[frm_idx + 1]
    session_dt_period = float(session.metadata(key='frm-dt-period'))
    session_tim = mbifile.datacube[f'frame-{frm_idx + 1}-data'][DataCubeKeys.AT_TIC][:]
    scan_entries = []
    num_gates = []
    num_samples = []
    unpacked = []
    scan_tic = []
    adc_scans = []
    gate_scans = []
    for scan_idx in range(len(session)):
        scan_tic.append(session_tim[scan_idx])
        scan = {}
        adc_idx = 0
        my_adc = session.data_counts[session.index_counts[scan_idx - 1]:
                                  session.index_counts[scan_idx]]
        my_gates = session.data_positions[session.index_positions[scan_idx - 1]:
                                       session.index_positions[scan_idx]]
        gates, adc = session.ms_raw(scan_idx + 1)
        num_gates.append(len(gates))
        num_samples.append(len(adc))
        #print('break')
        for gate in gates:
            for gate_idx in range(gate[0], gate[1]):
                scan[gate_idx] = adc[adc_idx]
                gate_scans.append(gate)
                adc_scans.append(adc[adc_idx])
                adc_idx += 1
                #print('break')
        if adc_idx != len(adc):
            raise RuntimeError('Scan gate and data size mismatch!')
        unpacked.append(scan)
    unpacked_dict.update({f'frame {frm_idx}': {'num_gates': num_gates, 'num_samples': num_samples,
                'scans': unpacked, 'scan_tic': scan_tic, 'timestamps': list(session.trigger_timestamps[:])}})
    at_gate_intensity_dict = {}
    for index, item in enumerate(unpacked):
        if item != {}:
            arrival_time = index*session_dt_period
            atgateintensity_entry = {arrival_time: item}
            at_gate_intensity_dict.update(atgateintensity_entry)
        else:
            pass

    temp_scan_holder = list(at_gate_intensity_dict.values())
    x_values = []
    y_values = []
    for thing in temp_scan_holder:
        x_values.extend(list(thing.keys()))
        y_values.extend(list(thing.values()))
    print('break')
    y = 0
    plt.xlabel('digitizer scans')
    plt.ylabel('intensity')
    for x in x_values:
        x_temp = [x, x]
        y_temp = [0, y_values[y]]
        y += 1
        plt.plot(x_temp, y_temp, c='blue')
    plt.show()
    print('break')

print('break')
#scan_entries = unpacked_dict.get('scans') #gate index values and intensity stored in a dictionary at a position
#num_gates = unpacked_dict.get('num_gates') # index of ToF events (as an array)
#num_samples = unpacked_dict.get('num_samples') # number of scans in each ToF event
#scan_tic = unpacked_dict.get('scan_tic') # time in ms(?) of each recording pulse