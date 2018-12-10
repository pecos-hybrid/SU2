inline void CAvgGrad_Hybrid::SetAniso_Eddy_Viscosity(su2double** aniso_eddy_viscosity_i,
                                                     su2double** aniso_eddy_viscosity_j) {
  Aniso_Eddy_Viscosity_i = aniso_eddy_viscosity_i;
  Aniso_Eddy_Viscosity_j = aniso_eddy_viscosity_j;
}


inline void CAvgGrad_Hybrid::SetPrimitive_Average(su2double* val_primvar_average_i,
                                                  su2double* val_primvar_average_j) {
  PrimVar_Average_i = val_primvar_average_i;
  PrimVar_Average_j = val_primvar_average_j;
}

/*!
 * \brief Set the gradient of the primitive variables for the average.
 * \param[in] val_primvar_grad_i - Gradient of the primitive variable at point i.
 * \param[in] val_primvar_grad_j - Gradient of the primitive variable at point j.
 */
inline void CAvgGrad_Hybrid::SetPrimVarGradient_Average(su2double **val_primvar_grad_i,
                                                        su2double **val_primvar_grad_j) {
  PrimVar_Grad_Average_i = val_primvar_grad_i;
  PrimVar_Grad_Average_j = val_primvar_grad_j;
}
