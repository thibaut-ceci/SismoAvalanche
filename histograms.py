def open_catalog(avalanches, pos_number=0.260, xlim=100):
    ax = avalanches.count().sort_values().plot(kind="barh", figsize=(6, 10))

    ax.set_xlabel("Number of avalanches")
    ax.set_xlim(0, xlim)

    #For counting the number of values in each column
    for index, value in enumerate(avalanches.count().sort_values()):
        ax.text(value+0.1, index-pos_number, str(value))


def see_word_distribution(ESEC, pos_number, xlim):
    ESEC_subtype = ESEC.value_counts().sort_values(ascending=False)

    ax = ESEC_subtype.plot(kind='barh')

    ax.set_xlabel("Number of avalanches")
    ax.set_yticks(range(len(ESEC_subtype)))
    ax.set_yticklabels(ESEC_subtype.index)
    ax.set_xlim(0, xlim)

    #For counting the number of values in each column
    for index, value in enumerate(ESEC.value_counts().sort_values(ascending=False)):
        ax.text(value+0.1, index-pos_number, str(value))


def see_number_distribution(avalanches_select, ax, ylabel):
    ax.hist(avalanches_select, orientation="horizontal")
    ax.set_xlabel("Number of avalanches")
    ax.set_ylabel(ylabel)
