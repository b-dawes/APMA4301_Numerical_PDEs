FILE(REMOVE_RECURSE
  "CMakeFiles/genwrappers"
  "SystemFunctionalsWrapper.cpp"
  "SystemSolversWrapper.cpp"
  "SystemExpressionsWrapper.cpp"
  "VisualizationWrapper.cpp"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/genwrappers.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
