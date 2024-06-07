import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd




#colors used in paper so far
refitting_colors_1 = ['#767B7B','#0F5B66', '#077E97', '#19874B', '#0ab8db',
                      '#AB242C','#54102E']

# green and blue scale
refitting_colors_2 = ["#0f3f12", "#328437", "#56a35b", "#6ddb74", "#90EE90", 
                      '#0066aa', '#57a6f3']

#rainbow scale with apobec signatures in purple and pink
refitting_colors_3 = ['#970707', '#974f07', '#979707', '#4f9707', '#077e97',
                      '#4f0797', '#97074f']

custom_colors_in = refitting_colors_2

#0F5B66 dark greenblue
#077E97 light green blue

def grouped_barplot_with_points(input_df):
    plt.figure(figsize=(4, 4.5))
    # Define custom colors
    custom_palette = {"SBS2": "#57a6f3", "SBS13": "#0066aa"}
    custom_palette_2 = {"SBS2": "#077E97", "SBS13": "#0F5B66"}
    # Draw the bar chart
    ax = sns.barplot(
        data=input_df, 
        x="condition", 
        y="value", 
        hue="signature", 
        alpha=0.7, 
        ci=None,
        palette=custom_palette_2  # Assign the custom palette
    )

    # Get the legend from just the bar chart
    handles, labels = ax.get_legend_handles_labels()

    # Draw the stripplot
    sns.stripplot(
        data=input_df, 
        x="condition", 
        y="value", 
        hue="signature", 
        dodge=True, 
        edgecolor="black", 
        linewidth=.75,
        ax=ax,
        palette=custom_palette_2  # Assign the custom palette
    )

    # Remove the old legend
    ax.legend_.remove()

    # Add just the bar chart legend back
    ax.legend(
        handles,
        labels,
    )
    # Add y-axis label
    ax.set_ylabel('Relative Contribution')

    return plt

#########################################################

def plot_grouped_replicates(input_data, plot_title:str, y_label:str):

    # Define custom colors for the legend
    custom_colors = custom_colors_in

    # Custom labels for x-axis
    x_labels = ['day 6', 'day 10', 'day 14',
                'day 17', 'day 23', 'day 28']

    #bar positions
    bar_positions = np.array([0.2, 0.6, 1,
                            0.2, 0.6, 1,
                            0.2, 0.6, 1,
                            0.2, 0.6, 1,
                            0.2, 0.6, 1,
                            0.2, 0.6, 1])

    # Plot the stacked bar plot
    ax = input_data.plot(x='sample', kind='bar',
                         stacked=True, color=custom_colors,
                         position=bar_positions)

    # Add y-axis label
    ax.set_ylabel(y_label)

    # Add title
    ax.set_title(plot_title)


    plt.legend(title = 'Signatures', loc='upper left',
               bbox_to_anchor=(1.012, 0.7))

    # Add vertical lines to separate groups of bars
    for i in range(1, len(x_labels)):
        plt.axvline(x=i*3 - 0.5, color='gray', linestyle='--')

    # Set x-tick positions and labels
    ax.set_xticks([1, 4, 7, 10, 13, 16])
    ax.set_xticklabels(x_labels, rotation=0)
    ax.set_xticklabels(x_labels)

    return plt
#####################################################
def plot_genotype_comparison_replicates(input_data, plot_title:str, y_label:str):

    # Define custom colors for the legend
    custom_colors = custom_colors_in

    # Custom labels for x-axis
    x_labels = ['control', 'sgUBR4',
                'sgUBR5', 'sgHUWE1']

    #bar positions
    bar_positions = np.array([0.2, 0.6, 1,
                            0.2, 0.6, 1,
                            0.2, 0.6, 1,
                            0.2, 0.6, 1])

    # Plot the stacked bar plot
    ax = input_data.plot(x='sample', kind='bar',
                         stacked=True, color=custom_colors,
                         position=bar_positions)

    # Add y-axis label
    ax.set_ylabel(y_label)

    # Add title
    ax.set_title(plot_title)


    plt.legend(title = 'Signatures', loc='upper left',
               bbox_to_anchor=(1.012, 0.7))

    # Add vertical lines to separate groups of bars
    for i in range(1, len(x_labels)):
        plt.axvline(x=i*3 - 0.5, color='gray', linestyle='--')

    # Set x-tick positions and labels
    ax.set_xticks([1, 4, 7, 10])
    ax.set_xticklabels(x_labels, rotation=0)
    ax.set_xticklabels(x_labels)

    return(plt)
####################################################
def plot_genotype_comparison_summary(input_data, plot_title:str):

    # Define custom colors for the legend
    custom_colors = custom_colors_in

    # Custom labels for x-axis
    x_labels = ['control', 'sgUBR4',
                'sgUBR5', 'sgHUWE1']
    
    #bar positions
    bar_positions = np.array([4, 3, 2, 1])

    # Plot the stacked bar plot
    ax = input_data.plot(x='sample', kind='bar',
                         stacked=True, color=custom_colors,
                         width = 5, position = bar_positions)

    # Add y-axis label
    ax.set_ylabel('Relative Contribution')

    # Add title
    ax.set_title(plot_title)


    plt.legend(title = 'Signatures', loc='upper left',
               bbox_to_anchor=(1.012, 0.7))

    # Set x-tick positions and labels
    ax.set_xticks([-18, -12, -6, 0])
    ax.set_xticklabels(x_labels, rotation=0)
    ax.set_xticklabels(x_labels)

    return plt

    #################################################
def normalize_refitting(df, condition: str):
    df_sub = df[df['condition'] == condition]
    df_sub = df_sub.iloc[:,0:7]
    #normalize to rowsums
    row_sums = df_sub.sum(axis=1)
    df_sub = df_sub.div(row_sums, axis= 0)
    # add and rename index
    df_sub = df_sub.reset_index()
    df_sub = df_sub.rename(columns={'index' : 'sample'})

    #change order of columns
    df_sub = df_sub[['sample', 'SBS1', 'SBS5', 'SBS6', 'SBS15', 'SBS18', 'SBS13', 'SBS2']]
    
    return df_sub
