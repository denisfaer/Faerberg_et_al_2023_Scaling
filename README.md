### Code repository for Faerberg, Gurarie & Ruvinsky *Dev Biol* 2026

Folders labeled **appendix2** and **figure5** provide codes to recapitulate analysis in the relevant figures.

The sample CSV sheet ``C_elegans_embryos.csv`` is a sample dataset for ``bootsrap.py`` and is a courtesy of Dr. Bao (Santella et al., 2010 *BMC Bioinform*; Moore et al., 2013 *Cell*; Du et al., 2014 *Cell*).

The Python script ``bootsrap.py`` is as an easy way to perform the CV analysis and get an empirical p-value.

### How to run the ``bootsrap.py`` script:

1. Format your data into a basic CSV file where rows represent individuals and columns represent timepoints. This file shouldn't contains labels or anything else besides raw data. See the ``C_elegans_embryos.csv`` file for reference.
2. Download the ``bootsrap.py`` code and place it in the same folder as your formatted CSV file.
3. Change the string in line 90 of ``bootsrap.py`` replacing *file_name_here* with the name of the formatted CSV file containing your data: make sure to keep the single quotation marks and DO NOT include the .csv extention.
4. Run the ``bootsrap.py`` code. It should generate the otput CSV files in the same directory as well as a new "Results" folder with the relevant bar plots.

If any issues arise try making a new environment and installing the following versions of Python and libraries:

``Python 3.10.4``
``numpy 1.21.5``
``scipy 1.7.3``
``matplotlib 3.5.3``
``seaborn 0.12.2``
