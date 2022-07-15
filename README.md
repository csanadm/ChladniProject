# ChladniProject
Fit data from Chladni plots with harmonic curves and plot them. For input the code needs a JSON file from such as the one produced with [Automeris WebPlotDigitizer](https://automeris.io/WebPlotDigitizer/). It furthermore requires gcc (tested with 5.4), GNU Make (tested with 4.1), ROOT (tested with 6.15), the latter also including the Minuit2 component (for the fitting).

## File content
- [**README.md**](https://github.com/csanadm/ChladniProject/blob/main/README.md): This README file
- [**Makefile**](https://github.com/csanadm/ChladniProject/blob/main/Makefile): Using `make <basename>.exe` or `make all`, it will create an executable from any `.cc` file
- [**json.hpp**](https://github.com/csanadm/ChladniProject/blob/main/json.hpp): JSON reader from [github.com/nlohmann/json](https://github.com/nlohmann/json/)
- [**chladni_fits.cc**](https://github.com/csanadm/ChladniProject/blob/main/chladni_fits.cc): Code for fitting and plotting
- [**shu2022_wpd_fig7.json**](https://github.com/csanadm/ChladniProject/blob/main/shu2022_wpd_fig7.json): Example data file from the 7th subplot of Fig. 2(b) of [Shu et al., Entropy 2022, 24(2), 215](https://www.mdpi.com/1099-4300/24/2/215), data extracted with [Automeris WebPlotDigitizer](https://automeris.io/WebPlotDigitizer/)
- [**shu2022_wpd_fig7_all_fit.png**](https://github.com/csanadm/ChladniProject/blob/main/shu2022_wpd_fig7_all_fit.png): Resulting plot if using the above example JSON file.

## Example output
![Example output plot]([https://github.com/adam-p/markdown-here/raw/master/src/common/images/icon48.png](https://raw.githubusercontent.com/csanadm/ChladniProject/main/shu2022_wpd_fig7_all_fit.png) "Example output plot")
