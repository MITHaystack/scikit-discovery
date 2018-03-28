class DerampSentinel(object):
    def __init__(self, metadata):

        self._metadata = metadata
        tree = self._metadata['Tree']

        burst_list = tree.findall('swathTiming/burstList/burst')
        az_fm_rate_list = tree.findall('generalAnnotation/azimuthFmRateList/azimuthFmRate')
        doppler_centroid_list = tree.findall('dopplerCentroid/dcEstimateList/dcEstimate')

        self._num_bursts = len(burst_list)


        self._vel_interp = OrbitInterpolation(self._metadata['Orbit'],interp_target='velocity')

        self._lines_per_burst = int(tree.find('swathTiming/linesPerBurst').text)
        self._samples_per_burst = int(tree.find('swathTiming/samplesPerBurst').text)

        self._az_time_interval = float(tree.find('imageAnnotation/imageInformation/azimuthTimeInterval').text)
        self._az_steering_rate = np.deg2rad(float(tree.find('generalAnnotation/productInformation/azimuthSteeringRate').text))

        radar_freq = float(tree.find('generalAnnotation/productInformation/radarFrequency').text)
        self._radar_lambda = c/radar_freq

        self._slant_range_time = float(tree.find('imageAnnotation/imageInformation/slantRangeTime').text)
        self._slant_range_time_interval = 1/float(tree.find('generalAnnotation/productInformation/rangeSamplingRate').text)

        self._doppler_centroid_scanning_rate_list = list(map(self._dopplerCentroidRate, burst_list))
        self._doppler_fm_rate_list = list(map(self._dopplerFMRate, az_fm_rate_list))
        self._doppler_centroid_frequency_list = list(map(self._dopplerCentroidFrequency, doppler_centroid_list))


#     def __call__(self, lines, samples, invert=False):

#         piecewise_functions = [lambda lines, samples: self._derampedPhase(lines, samples, index) for index in range(self._num_bursts)]

#         np.piecewise()


#    def _derampedPhase(self, lines, samples, index):
    def __call__(self, lines, samples, index):
        centroid_tops = self._dopplerCentroidTops(samples,
                                                  self._doppler_fm_rate_list[index],
                                                  self._doppler_centroid_scanning_rate_list[index])

        zero_doppler_azimuth_time = self._zeroDopplerAzimuthTime(lines)
        ref_zero_doppler_azimuth_time = self._referenceZeroDopplerAzimuthTime(samples,
                                                                              self._doppler_centroid_frequency_list[index],
                                                                              self._doppler_fm_rate_list[index])


        return -np.pi * centroid_tops * (zero_doppler_azimuth_time - ref_zero_doppler_azimuth_time)**2


    def _dopplerCentroidRate(self, burst):

        az_start_time = pd.to_datetime(burst.find('azimuthTime').text)

        az_time_mid_burst =   az_start_time \
                            + pd.to_timedelta(self._az_time_interval*self._lines_per_burst/2,'s')

        speed = np.linalg.norm(self._vel_interp(az_time_mid_burst))

        return self._az_steering_rate * 2 * speed / self._radar_lambda

    def _dopplerFMRate(self, az_fm_rate):
        doppler_fm_rate_t0 = float(az_fm_rate.find('t0').text)
        doppler_fm_rate_coeffs = [float(az_fm_rate.find(label).text) for label in ['c0','c1','c2']]

        return DerampPolynomial(doppler_fm_rate_t0,
                                doppler_fm_rate_coeffs,
                                self._slant_range_time_interval,
                                self._slant_range_time)

    def _dopplerCentroidTops(self, samples, doppler_fm_rate, doppler_centroid_rate_from_scanning_antenna):

        return (doppler_fm_rate(samples) * doppler_centroid_rate_from_scanning_antenna) \
              / (doppler_fm_rate(samples) - doppler_centroid_rate_from_scanning_antenna)



    def _dopplerCentroidFrequency(self, doppler_centroid_estimate):
        doppler_centroid_coeffs = [float(poly) for poly in doppler_centroid_estimate.find('dataDcPolynomial').text.split(' ')]
        doppler_centroid_t0 = float(doppler_centroid_estimate.find('t0').text)
        return DerampPolynomial(doppler_centroid_t0,
                          doppler_centroid_coeffs,
                          self._slant_range_time_interval,
                          self._slant_range_time)


    def _zeroDopplerAzimuthTime(self, lines):
        lines = lines % self._lines_per_burst

        return lines*self._az_time_interval - self._az_time_interval * self._lines_per_burst/2

    def _referenceZeroDopplerAzimuthTime(self, samples, doppler_centroid_frequency, doppler_fm_rate):
        beam_center_crossing_time = lambda samples: -doppler_centroid_frequency(samples) / doppler_fm_rate(samples)
        return beam_center_crossing_time(samples) - beam_center_crossing_time(0)
