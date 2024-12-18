set(SRC
  # amr
  amr/AMR.cpp
  amr/CurrentCalculator.cpp
  amr/LinearSolver.cpp

  # anisotropy
  anisotropy/anisotropy.cpp
  anisotropy/anisotropy_cpu.cpp
  anisotropy/fs_anisotropy.cpp
  anisotropy/fs_anisotropy_cpu.cpp
  anisotropy/anisotropy_new.cpp
  anisotropy/anisotropy_new_cpu.cpp


  #dmi
  dmi/dmi.cpp
  dmi/dmi_cpu.cpp
  dmi/fs_dmi.cpp
  dmi/fs_dmi_cpu.cpp

  # demag
  demag/phi/demag_phi_tensor.cpp
  demag/tensor_round.cpp
  demag/demag_static.cpp
  demag/tensor.cpp
  demag/old/demag_old.cpp
  demag/demag_tensor.cpp
  demag/fs_tensor.cpp
  demag/fs_demag_tensor.cpp
  demag/multiscale.cpp

  #temperature
  temperature/fdm_temperature.cpp
  temperature/fdm_temperature_cpu.cpp
  temperature/fs_fdm_temperature.cpp
  temperature/fs_fdm_temperature_cpu.cpp

  # exchange
  exchange/exchange.cpp
  exchange/exchange_cpu.cpp
  exchange/fs_exchange.cpp
  exchange/fs_exchange_cpu.cpp

  # exchange_nnn
  exchange_nnn/exchange_nnn.cpp
  exchange_nnn/exchange_nnn_cpu.cpp

  #spinhall
  spinhall/fs_spinhall.cpp
  spinhall/fs_spinhall_cpu.cpp

  # io
  io/OMFHeader.cpp
  io/OMFExport.cpp
  io/OMFImport.cpp

  # llge
  llge/llge.cpp
  llge/llge_cpu.cpp

  # spintorque
  spintorque/fdm_zhangli.cpp
  spintorque/fdm_zhangli_cpu.cpp

  # macrospintorque
  macrospintorque/fdm_slonchewski.cpp
  macrospintorque/fdm_slonchewski_cpu.cpp
  macrospintorque/fs_fdm_slonchewski.cpp
  macrospintorque/fs_fdm_slonchewski_cpu.cpp
  macrospintorque/fdm_slonchewski_angleDep.cpp
  macrospintorque/fdm_slonchewski_cpu_angleDep.cpp

  # charge
  topology/charge.cpp

  # minimize
  minimize/minimize.cpp
  minimize/minimize_cpu.cpp


  # interlayerExchange
  interlayerExchange/interlayerExchange.cpp
  interlayerExchange/interlayerExchange_cpu.cpp
  interlayerExchange/interlayerExchange_multi.cpp
  interlayerExchange/interlayerExchange_multi_cpu.cpp




)

set(CUSRC
  anisotropy/anisotropy_cuda.cu
  anisotropy/anisotropy_new_cuda.cu
  temperature/fdm_temperature_cuda.cu
  temperature/fs_fdm_temperature_cuda.cu
  exchange/exchange_cuda.cu
  exchange_nnn/exchange_nnn_cuda.cu
  llge/llge_cuda.cu
  spintorque/fdm_zhangli_cuda.cu
  anisotropy/fs_anisotropy_cuda.cu
  dmi/dmi_cuda.cu
  interlayerExchange/interlayerExchange_cuda.cu
  interlayerExchange/interlayerExchange_multi_cuda.cu
  dmi/fs_dmi_cuda.cu
  exchange/fs_exchange_cuda.cu
  macrospintorque/fs_fdm_slonchewski_cuda.cu
  macrospintorque/fdm_slonchewski_cuda.cu
  macrospintorque/fdm_slonchewski_cuda_angleDep.cu
  minimize/minimize_cuda.cu
)

append_magneto_src(${CMAKE_CURRENT_SOURCE_DIR} "${SRC}")
append_magneto_cusrc(${CMAKE_CURRENT_SOURCE_DIR} "${CUSRC}")
