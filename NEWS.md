# Version 0.0.1

### Initial release

- Initial release to CRAN.
- Changes to lablaster up to this point are documented at https://github.com/alexsb1/lablaster.

# Version: 1.0.1

### Major update during peer-review
- An argument used in the endPoint() function was changed from "df" to "detectDf" to avoid conflict with already published packages. This is a breaking change only if the argument name "df" is explicitly referenced. The order of arguments remain unchanged so this is a potential breaking change.
- Added over-smoothing warning.
- Made scanRate dynamic by fixing a grep of time column.
- The returned data frame now contains a data frame subset of the original data frame supplied that contains only observations less than endTime.
