import datetime
import os
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import shap

from core_xfp import FingerprintManager
from report_maker import BitAnalysis


class XFPTreeExplainer(FingerprintManager):
    def __init__(self, in_tree_explainer: shap.TreeExplainer) -> None:
        """
        Initialize the XFPTreeExplainer class.

        Parameters
        ----------
        in_tree_explainer : shap.TreeExplainer
            Tree explainer object from SHAP to use.

        Raises
        ------
        ValueError
            If the input tree explainer is not a valid SHAP Tree Explainer.
        """
        super().__init__()
        self.shap_values = None
        self.tree_explainer = in_tree_explainer

        # Check whether the tree explainer is valid or not
        if not isinstance(self.tree_explainer, shap.TreeExplainer):
            raise ValueError("Invalid SHAP Tree Explainer")

    def generate_shap_values(self, return_shap_values: bool = False) -> Optional[np.ndarray]:
        """
        Generate SHAP values of the input fingerprint using the Tree Explainer.

        Parameters
        ----------
        return_shap_values : bool, optional
            Return the SHAP values or not (default is False).

        Returns
        -------
        np.ndarray or None
            SHAP values of the input fingerprint in `self.fp` if `return_shap_values` is True, None otherwise.
        """
        self.shap_values = self.tree_explainer.shap_values(self.fp)

        if return_shap_values:
            return self.shap_values

    def _folder_settings(self) -> None:
        """
        Asks for input for folder (with path) to store the images and the report.
        Validates the path and creates the folder if it doesn't exist.
        """

        # Ask for input for folder (with path) to store the images and the report
        self.folder_name = input(
            f"Enter the folder name (with path) to store the images and the report. The default is 'XFP_Tree_Explainer_(today's date)' "
        )

        # Validate the path and create the folder if it doesn't exist
        if self.folder_name == "":
            # Add date and a folder name
            self.folder_name = f"XFP_Tree_Explainer_{datetime.date.today()}"

        if not os.path.exists(self.folder_name):
            os.makedirs(self.folder_name)

    def generate_shap_summary_plot(self, max_display: int = 10, save_fig: str = "", dpi: int = 300) -> None:
        """
        Generate a SHAP summary plot.

        Parameters
        ----------
        max_display : int, optional
            Maximum number of features to plot (default is 10).
        save_fig : str, optional
            If provided, save the figure instead of showing it (default is "").
        dpi : int, optional
            Dots per inch (DPI) for the figure (default is 300).
        """

        # Feature names
        feature_names = [f"Bit {i}" for i in range(self.fp.shape[1])]

        # Check if save_fig is provided, if yes, save the figure instead of showing it.
        if save_fig != "":
            shap.summary_plot(
                self.shap_values,
                self.fp,
                max_display=max_display,
                feature_names=feature_names,
                show=False,
            )
            plt.savefig(save_fig, dpi=dpi, bbox_inches="tight")
            plt.close()
            return

        shap.summary_plot(self.shap_values, self.fp, max_display=max_display, feature_names=feature_names)

    def generate_single_shap_summary_plot(self, query_bit: int, return_fig_path: bool = True) -> Optional[str]:
        """
        Saves a single SHAP summary plot of the SHAP values of the query bit.

        Parameters
        ----------
        query_bit : int
            Bit for which the SHAP Summary Plot is to be generated.
        return_fig_path : bool, optional
            Return the path of the figure or not (default is True).

        Returns
        -------
        str or None
            Path of the figure if `return_fig_path` is True, None otherwise.
        """

        # Making the summary plot of only the query bit.
        shap.summary_plot(
            self.shap_values[:, query_bit : query_bit + 1],
            self.fp[:, query_bit : query_bit + 1],
            feature_names=[f"Bit {query_bit}"],
            show=False,
        )

        # Save the plot
        shap_summary_plot_path = f"{self.folder_name}/Bit_{query_bit}_SHAP_Summary_Plot.png"
        plt.savefig(shap_summary_plot_path, dpi=1000, bbox_inches="tight")

        plt.clf()
        plt.close()

        if return_fig_path:
            return shap_summary_plot_path

    def generate_bit_analysis_report(
        self,
        bit_list: list[int],
        report_title: str = "",
        time_stamp: bool = True,
        file_name: str = "X-FP_Bit_Analysis_Report",
    ) -> None:
        """
        Generate a PDF report for the analysis of the query bits.

        Parameters:
        -----------
        bit_list : list[int]
            List of bits for which the analysis is to be performed.
        report_title : str, optional
            Title of the report provided by the user (default is "").
        time_stamp : bool, optional
            Add time stamp to the report footer or not (default is True). The time is local and is formatted as
            '%Y-%m-%d %H:%M:%S'.
        file_name : str, optional
            Name of the report file (default is "X-FP_Bit_Analysis_Report").
        """

        # folder settings
        self._folder_settings()

        # Start the report
        bit_analysis_report = BitAnalysis()

        # Generate the overall feature importance plot
        self.generate_shap_summary_plot(save_fig=f"{self.folder_name}/Overall_Feature_Importance_Plot.png")
        overall_shap_summary_plot_path = f"{self.folder_name}/Overall_Feature_Importance_Plot.png"

        # Print the title page
        bit_analysis_report.intro_page(
            overall_feature_importance_image=overall_shap_summary_plot_path,
            img_text="SHAP Summary Plot for the Top 10 Morgan Fingerprint Bits",
            report_title=report_title,
            time_stamp=time_stamp,
        )

        # Looping over all the query bits
        for query_bit in bit_list:
            # Generate the SHAP summary plot for the query bit
            shap_summary_plot_path = self.generate_single_shap_summary_plot(query_bit, return_fig_path=True)

            # Generating the count of substructures for the query bit, and retrieving the dataframe
            substructure_count_df = self.substruct_counter_per_bit(query_bit)

            # Generating the substructure images for the query bit, and saving that image in the folder
            substructure_images = self.get_substruct_images_per_bit(query_bit)
            substructure_images_path = f"{self.folder_name}/Bit_{query_bit}_Substructure_Images.svg"
            with open(substructure_images_path, "w") as svg_file:
                svg_file.write(substructure_images)

            # Generating the bit analysis report
            bit_analysis_report.print_chapter(
                in_bit=query_bit,
                mol_importance_image=shap_summary_plot_path,
                mol_image=substructure_images_path,
                in_df=substructure_count_df,
                n_compounds=len(self.input_df),
            )

        # Add end note to the report
        bit_analysis_report.add_end_note(feature_importance_note="shap_tree_explainer")

        # Save the report
        bit_analysis_report.output(f"{self.folder_name}/{file_name}.pdf")
