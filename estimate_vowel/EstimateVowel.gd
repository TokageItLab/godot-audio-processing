extends Control

enum Vowel {
    A,
    I,
    U,
    E,
    O
}

const UPDATE_FRAME: int = 5 # getting the recording stream again at a rate of 4 frames or less caused error, WTF
const FFT_SAMPLES: int = 1024 # needs 2^n

const GRAPH_LENGTH: float = 640.0
const GRAPH_SCALE: float = 200.0

const DYNAMIC_RANGE: float = 100.0 # assuming the min value of db
const INV_DYNAMIC_RANGE: float = 1.0 / DYNAMIC_RANGE # for reducing calculation
const PI2: float = 2.0 * PI
const INV_255: float = 1.0 / 255.0
const INV_32767: float = 1.0 / 32767.0
const INV_LOG10: float = 1.0 / log(10)

# Vector2(freq, value ratio)
const ESTIMATE_DB: Dictionary = {
    "peak3": {
        "A": [
            Vector2(18, 1.0),
            Vector2(41, 0.9),
            Vector2(85, 0.75)
        ],
        "I": [
            Vector2(21, 1.0),
            Vector2(42, 1.1),
            Vector2(84, 1.0)
        ],
        "U": [
            Vector2(19, 1.0),
            Vector2(47, 0.65),
            Vector2(84, 0.7)
        ],
        "E": [
            Vector2(21, 1.0),
            Vector2(60, 0.75),
            Vector2(84, 0.65)
        ],
        "O": [
            Vector2(20, 1.0),
            Vector2(63, 0.9),
            Vector2(85, 0.8)
        ]
    },
    "peak4": {
        "A": [
            Vector2(18, 1.0),
            Vector2(41, 0.9),
            Vector2(68, 0.7),
            Vector2(85, 0.55)
        ],
        "I": [
            Vector2(21, 1.0),
            Vector2(42, 1.1),
            Vector2(60, 1.0),
            Vector2(84, 1.1)
        ],
        "U": [
            Vector2(20, 1.0),
            Vector2(39, 0.7),
            Vector2(65, 0.6),
            Vector2(84, 0.75)
        ],
        "E": [
            Vector2(22, 1.0),
            Vector2(43, 0.9),
            Vector2(66, 0.7),
            Vector2(84, 0.65)
        ],
        "O": [
            Vector2(20, 1.0),
            Vector2(39, 0.9),
            Vector2(63, 0.75),
            Vector2(85, 0.8)
        ]
    }
}

var effect: AudioEffectRecord
var audio_sample: AudioStreamSample
var is_recording: bool = false
var buffer: int = UPDATE_FRAME

var before_sample_array: Array = []
var peaks3_log: Array = []
var peaks4_log: Array = []
var vowel_log: Array = [-1, -1, -1]
var estimate_log: Array = [-1, -1, -1]

#########
# debug #
#########
# get and print max value in array
func print_max(sample_array: Array):
    var vmax: float = 0.0;
    for i in range(sample_array.size()):
        vmax = max(vmax, abs(sample_array[i]))
    print(vmax)
    return
func print_min(sample_array: Array):
    var vmin: float = 1.0;
    for i in range(sample_array.size()):
        vmin = min(vmin, abs(sample_array[i]))
    print(vmin)
    return

# (Caution!) plotting many points cause unstable.
func draw_graph(sample_array: Array) -> void:
    $Line2D.clear_points()
    var step: int = 0;
    for i in range(sample_array.size()):
        var px = (float(GRAPH_LENGTH) / float(sample_array.size())) * float(i)
        var py = float(-sample_array[i]) * float(GRAPH_SCALE)
        $Line2D.add_point(Vector2(px, py))
        step += 1
    return

# text = vowel, opacity = amount 
func draw_vowel(vowel: int, amount: float) -> void:
    var text: String = ""
    match vowel:
        Vowel.A:
            text = "A"
        Vowel.I:
            text = "I"
        Vowel.U:
            text = "U"
        Vowel.E:
            text = "E"
        Vowel.O:
            text = "O"
        _:
            text = "-"
    
    $Label.text = text
    $Label.set_modulate(Color(1,1,1,amount))
    return

#####################
# general functions #
#####################

# read mic input sample
# reference (https://godotengine.org/qa/67091/how-to-read-audio-samples-as-1-1-floats) 
static func read_16bit_samples(stream: AudioStreamSample) -> Array:
    assert(stream.format == AudioStreamSample.FORMAT_16_BITS)
    var bytes = stream.data
    var samples = []
    var i = 0
    # Read by packs of 2 bytes
    while i < len(bytes):
        var b0 = bytes[i]
        var b1 = bytes[i + 1]
        # Combine low bits and high bits to obtain 16-bit value
        var u = b0 | (b1 << 8)
        # Emulate signed to unsigned 16-bit conversion
        u = (u + 32768) & 0xffff
        # Convert to -1..1 range
        var s = float(u - 32768) / 32768.0
        samples.append(s)
        i += 2
    return samples

# get rms
func calc_rms(sample_array: Array) -> float:
    var rms: float = 0.0
    for i in range(sample_array.size()):
        rms += sample_array[i] * sample_array[i]
    rms = sqrt(rms / sample_array.size())
    rms = 20 * (log(rms) * INV_LOG10)
    return rms

# normalize
func array_normalize(sample_array: Array):
    var n: int = sample_array.size()
    var vmax: float = 0.0;
    var vmin: float = 0.0;
    for i in range(sample_array.size()):
        vmax = max(vmax, sample_array[i])
        vmin = min(vmin, sample_array[i])
    var d: float = 1.0 / (vmax - vmin) if (vmax - vmin) != 0 else 1.0
    for i in range(n):
        sample_array[i] = (sample_array[i] - vmin) * d
    return

#####################
# functions for FFT #
#####################

# smoothing
func smoothing(sample_array: Array):
    var n = sample_array.size();
    for i in range(n):
        sample_array[i] = (sample_array[i] + before_sample_array[i]) * 0.5
    return

# hamming window
func hamming(sample_array: Array):
    var n = sample_array.size();
    for i in range(n):
        var h = 0.54 - 0.46 * cos(PI2 * i / float(n - 1));
        sample_array[i] = sample_array[i] * h;
    sample_array[0] = 0
    sample_array[n - 1] = 0
    return;

# real FFT (wrap fft())
func rfft(sample_array: Array, reverse: bool = false, positive: bool = true):
    var n: int = sample_array.size()
    var cmp_array = []
    for i in range(n):
        cmp_array.append(Vector2(sample_array[i], 0.0))
    fft(cmp_array, reverse)
    if positive:
        for i in range(n):
            sample_array[i] = abs(cmp_array[i].x)
    else:
        for i in range(n):
            sample_array[i] = cmp_array[i].x
    if reverse:
        var inv_n: float = 1.0 / float(n)
        for i in range(n):
            sample_array[i] *= inv_n
    return

# FFT
# reference (https://caddi.tech/archives/836) 
func fft(a: Array, reverse: bool):
    var N: int = a.size()
    if N == 1: return
    var b: Array = []
    var c: Array = []
    for i in range(N):
        if i % 2 == 0:
            b.append(a[i])
        if i % 2 == 1:
            c.append(a[i])
    fft(b, reverse);
    fft(c, reverse);
    var circle: float = -PI2 if reverse else PI2
    for i in range(N):
        a[i] = b[i % (N / 2)] + ComplexCalc.cmlt(c[i % (N / 2)], ComplexCalc.cexp(Vector2(0, circle * float(i) / float(N))));
    return

# lifter
func lifter(sample_array: Array, level: int):
    var i_min: int = level
    var i_max: int = sample_array.size() - 1 - level
    for i in range(sample_array.size()):
        if i > i_min && i <= i_max:
            sample_array[i] = 0.0
    return

# filter
func filter(sample_array: Array, lowcut: int, highcut: int):
    var minimum = sample_array[0]
    for i in range(sample_array.size()):
        minimum = min(minimum, sample_array[i])
    if minimum == 0.0:
        minimum == 0.000001 # avoid log(0)
    for i in range(sample_array.size()):
        if sample_array[i] <= lowcut || sample_array[i] >= highcut:
            sample_array[i] = minimum

##########################
# functions for Lip Sync #
##########################

# get peaks
func get_peaks(sample_array: Array, threshold: float) -> Array:
    var n: int = sample_array.size() - 1
    var i: int = 1
    var j: int = 0
    var tmp: Vector2 = Vector2.ZERO
    var out: Array = []
    var div: float = 1.0
    while i < n:
        if (
            (sample_array[i] > threshold) &&
            (sample_array[i] > sample_array[i - 1]) &&
            (sample_array[i] > sample_array[i + 1])
        ):
            if out.size() > 0:
                out.append(Vector2(i, sample_array[i] * div))
            else:
                out.append(Vector2(i, 1.0))
                div = 1.0 / sample_array[i]
        i += 1
    return out

# get peaks average
func get_peaks_average(size: int) -> Array:
    var out: Array = []
    var i: int = 1
    var j: int = 0
    var div: float = 1.0
    match size:
        3:
            out = peaks3_log[0]
            while i < peaks3_log.size():
                j = 0
                while j < out.size():
                    out[j].x += peaks3_log[i][j].x
                    out[j].y += peaks3_log[i][j].y
                    j += 1
                i += 1
            div = 1.0 / peaks3_log.size()
        4:
            out = peaks4_log[0]
            while i < peaks4_log.size():
                j = 0
                while j < out.size():
                    out[j].x += peaks4_log[i][j].x
                    out[j].y += peaks4_log[i][j].y
                    j += 1
                i += 1
            div = 1.0 / peaks4_log.size()
    for k in range(out.size()):
        out[k] *= div
    return out

# get distance from vowel DB
func get_distance_from_db(peaks: Array) -> Array:
    var out: Array = []
    var vowel: Array = ["A", "I", "U", "E", "O"]
    var peak_estm: Dictionary = {}
    var dist: float = 0.0
    match(peaks.size()):
        3:
            peak_estm = ESTIMATE_DB["peak3"]
        4:
            peak_estm = ESTIMATE_DB["peak4"]
    for i in range(vowel.size()):
        dist = 0.0
        for j in range(peaks.size()):
            dist += abs(peak_estm[vowel[i]][j].x - peaks[j].x) * INV_255 + abs(peak_estm[vowel[i]][j].y - peaks[j].y)
        out.append(dist)
    return out

# estimate vowel by peak (too roughly)
func estimate_vowel(sample_array: Array) -> int:
    var out: Array = []
    var peaks: Array = get_peaks(sample_array, 0.1)

    if peaks.size() != 3 && peaks.size() != 4:
        return -1    
    push_peaks(peaks)
    var peaks_ave: Array = get_peaks_average(peaks.size())
    print(peaks_ave)
    
    var distance_vowel: Array = get_distance_from_db(peaks_ave)
    print(distance_vowel)

    var i: int = 1
    var min_distance: float = distance_vowel[0]
    var min_idx = 0
    while i < 5:
        if distance_vowel[i] < min_distance:
            min_distance = distance_vowel[i]
            min_idx = i
        i += 1

    # unknown
    return min_idx

# estimate and complement (wrap estimate_vowel()
func get_vowel(sample_array: Array, amount: float) -> Dictionary:
    var current: int = estimate_vowel(sample_array)

    # input begin and end
    if vowel_log[0] != -1:
        if amount < 0.5:
            return {"estimate": current, "vowel": vowel_log[0] if vowel_log[0] != -1 else randi() % 5}

    # stabilize
    # now if (current == estimate_log[0]) will stabilize while 1 frame
    # so if (current == estimate_log[0] && current == estimate_log[1]) will stabilize while 2 frame
    if vowel_log.size() > 2:
        if current == estimate_log[0]:
            return {"estimate": current, "vowel": current if current != -1 else randi() % 5}
        else:
            return {"estimate": current, "vowel": vowel_log[0] if vowel_log[0] != -1 else randi() % 5}
    
    return {"estimate": current, "vowel": randi() % 5}

# log
func push_peaks(peaks: Array) -> void:
    if peaks.size() >= 4:
        if peaks4_log.size() < 3:
            peaks4_log.append(peaks)
        else:
            peaks4_log.push_front(peaks)
            peaks4_log.pop_back()
    elif peaks.size() >= 3:
        if peaks3_log.size() < 3:
            peaks3_log.append(peaks)
        else:
            peaks3_log.push_front(peaks)
            peaks3_log.pop_back()
    return
func push_vowel(vowel: int) -> void:
    if vowel_log.size() < 3:
        vowel_log.append(vowel)
    else:
        vowel_log.push_front(vowel)
        vowel_log.pop_back()
    return
func push_estimate(vowel: int) -> void:
    if estimate_log.size() < 3:
        estimate_log.append(vowel)
    else:
        estimate_log.push_front(vowel)
        estimate_log.pop_back()
    return

###########
# process #
###########

func _ready():
    var idx: int = AudioServer.get_bus_index("Record")
    effect = AudioServer.get_bus_effect(idx, 0)

func _on_RecordButton_pressed():
    is_recording = !is_recording
    if is_recording == false:
        $RecordButton.text = "Record"
    else:
        $RecordButton.text = "Stop"

func _process(delta):
    if is_recording:
        if buffer <= 0:
            if effect.is_recording_active():
                effect.set_recording_active(false)
                audio_sample = effect.get_recording()
                audio_sample.set_format(AudioStreamSample.FORMAT_16_BITS)
                if audio_sample:
                    var src_array: PoolByteArray = effect.get_recording().get_data()
                    var sample_array: Array = read_16bit_samples(audio_sample)
                    # calc RMS
                    var rms: float = calc_rms(sample_array)
                    # FFT
                    if sample_array.size() >= FFT_SAMPLES:
                        sample_array = sample_array.slice(0, FFT_SAMPLES - 1)
                        # hamming
                        hamming(sample_array)
                        # to spectrum by FFT
                        rfft(sample_array)
                        sample_array = sample_array.slice(0, FFT_SAMPLES * 0.5)
                        # smoothing
                        if !before_sample_array.empty():
                            smoothing(sample_array)
                        before_sample_array = sample_array.duplicate()
                        # filter
                        filter(sample_array, 10, 95) # adjust for formant count
                        # log power scale
                        for i in range(sample_array.size()):
                            sample_array[i] = log(pow(sample_array[i], 2)) * INV_LOG10
                        array_normalize(sample_array)
                        # to cepstrum by IFFT
                        rfft(sample_array, true, false)
                        # lifter
                        lifter(sample_array, 26) # adjust for formant count
                        # to spectrum by FFT again
                        rfft(sample_array, false, false)
                        sample_array = sample_array.slice(0, FFT_SAMPLES * 0.25)
                        # normalize
                        array_normalize(sample_array)
                        # emphasis peak
                        for i in range(sample_array.size()):
                            sample_array[i] = pow(sample_array[i], 2)
                        # normalize again and multiply RMS
                        array_normalize(sample_array)
                        var nrm_rms = min(DYNAMIC_RANGE, max(rms + DYNAMIC_RANGE, 0))
                        for i in range(sample_array.size()):
                            sample_array[i] = (sample_array[i] * nrm_rms * INV_DYNAMIC_RANGE)
                        # estimate vowel
                        var amount = clamp(inverse_lerp(-DYNAMIC_RANGE, 0.0, rms), 0.0, 1.0)
                        var current_vowel: Dictionary = get_vowel(sample_array, amount)
                        push_estimate(current_vowel["estimate"])
                        push_vowel(current_vowel["vowel"])
                        # draw
                        draw_vowel(current_vowel["vowel"], amount)
                        draw_graph(sample_array)
            effect.set_format(AudioStreamSample.FORMAT_16_BITS)
            effect.set_recording_active(true)
            buffer = UPDATE_FRAME
        else:
            buffer -= 1
    else:
        effect.set_recording_active(false)
    return
