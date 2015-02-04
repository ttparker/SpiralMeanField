std::vector<obsMatrixD_t, Eigen::aligned_allocator<obsMatrixD_t>> obsList(3);
obsList[0] << 0., 1.,
              1., 0.;
obsList[1] << std::complex<double>(0., 0.), std::complex<double>(0., -1.),
              std::complex<double>(0., 1.), std::complex<double>(0., 0.);
obsList[2] << 1., 0.,
              0., -1.;
