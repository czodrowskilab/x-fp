import datetime
import importlib.resources as pkg_resources
import os
from typing import Optional

import matplotlib.pyplot as plt
import pandas as pd
import xgboost as xgb
from xgboost import XGBClassifier

from . import resources
from .core_xfp import FingerprintManager
from .report_maker import BitAnalysis


class XGBFP(FingerprintManager):
    def __init__(self, in_xgboost_model: xgb.XGBClassifier) -> None:
        """
        Initialize the XGBFP class.

        Parameters
        ----------
        in_xgboost_model : xgboost.XGBClassifier
            XGBoost Classifier model to use.

        Raises
        ------
        ValueError
            If the input model is not a valid XGBoost Classifier model.
        """
        super().__init__()
        self.xgboost_model = in_xgboost_model

        # Check whether the model is a valid XGBoost Classifier model or not
        if not isinstance(self.xgboost_model, xgb.XGBClassifier):
            raise ValueError("Invalid object. Please input valid XGBoost Classifier model")

    def _folder_settings(self) -> None:
        """
        Asks for input for folder (with path) to store the images and the report.
        Validates the path and creates the folder if it doesn't exist.
        """

        # Ask for input for folder (with path) to store the images and the report
        self.folder_name = input(
            f"Enter the folder name (with path) to store the images and the report. The default is 'XFP_XGBoost_Classifier_(today's date)' "
        )

        # Validate the path and create the folder if it doesn't exist
        if self.folder_name == "":
            # Add date and a folder name
            self.folder_name = f"XFP_XGBoost_Classifier_{datetime.date.today()}"

        if not os.path.exists(self.folder_name):
            os.makedirs(self.folder_name)

    @staticmethod
    def _get_feature_importance(
        in_model: XGBClassifier,
        importance_type: str = "weight",
        n_bits: int = 2048,
    ) -> list[float]:
        """
        Get feature importance of the XGBoost Classifier model.

        Parameters
        ----------
        in_model : xgboost.XGBClassifier
            XGBoost Classifier model
        importance_type : str, optional
            Type of importance to be calculated (default is "weight"). The options are "weight", "gain", "cover",
            "total_gain", and "total_cover".
        n_bits : int, optional
            Number of bits in the fingerprint (default is 2048).

        Returns
        -------
        list[float]
            List of feature importance values for all the bits in the fingerprint.

        Notes
        -----
        This function is necessary because there isn't a straight-forward way to get feature importance from a XGBoost
        model for all (query) features. XGBoost doesn't store zero-importance features, plus there are different feature
        importance types to be considered.
        """

        # Assign names to the bits in the model object
        bit_names = [f"Bit {bit}" for bit in range(n_bits)]
        in_model.get_booster().feature_names = bit_names

        # Get feature importance
        importance_score = []
        for bit_name in bit_names:
            try:
                # XGBoost doesn't store zero-importance features
                importance_score.append(in_model.get_booster().get_score(importance_type=importance_type)[bit_name])
            except KeyError:
                importance_score.append(0)

        return importance_score

    def make_feature_importance_df(
        self,
        importance_type: str = None,
        calculate_all_importances: bool = False,
    ) -> None:
        """
        Make a DataFrame containing queried XGBoost feature importance.

        Parameters
        ----------
        importance_type : str, optional
            Type of importance to be calculated (default is None). The options are "weight", "gain", "cover",
            "total_gain", and "total_cover".
        calculate_all_importances : bool, optional
            Additionally calculate all the importances (default is False).

        Raises
        ------
        ValueError
            If the user doesn't provide either importance_type or calculate_all_importances or if the provided
            importance_type is invalid.
        """
        self.importance_df = pd.DataFrame()

        # validate if the user entered either importance type or calculate_all_importances
        if importance_type is None and not calculate_all_importances:
            raise ValueError("Please either enter an importance type or calculate all the importances.")

        # validate if the user entered correct importance type
        if importance_type not in ["weight", "gain", "cover", "total_gain", "total_cover"]:
            raise ValueError("Invalid importance type. Please enter a valid importance type.")

        self.importance_df["Bits"] = [f"Bit {bit}" for bit in range(self.n_bits)]
        if calculate_all_importances:
            importance_types = ["weight", "gain", "cover", "total_gain", "total_cover"]
            print(f"Calculating all importances: {importance_types}. Please wait.")
            for importance_type in importance_types:
                self.importance_df[importance_type] = self._get_feature_importance(
                    self.xgboost_model, importance_type=importance_type, n_bits=self.n_bits
                )
            print("Done.")

        else:
            self.importance_df[importance_type] = self._get_feature_importance(
                self.xgboost_model, importance_type=importance_type, n_bits=self.n_bits
            )

        if importance_type:
            self.primary_importance_type = importance_type

    def generate_feaure_importance_plot(
        self, importance_type: str = None, max_importance: int = 10, save_fig: str = ""
    ) -> None:
        """
        Generate the feature importance plot.

        Parameters
        ----------
        importance_type : str, optional
            Type of importance to be used for plotting (default is None).
            If the primary_importance_type is already set, then that will be used for plotting.
        max_importance : int, optional
            Maximum number of bits to be plotted (default is 10).
        save_fig : str, optional
            Path to save the figure (default is "", i.e. the image is not saved).

        Raises
        ------
        ValueError
            If the user hasn't run the `.make_feature_importance_df()` method first, or if the user doesn't provide an
            importance type or has not set the primary_importance_type, or if the user provides an invalid importance
            type.

        Notes
        -----
        Displays the plot if save_fig is not provided.
        """

        # validate if self.importance_df exists
        if self.importance_df is None:
            raise ValueError("Please run the make_feature_importance_df() method first.")

        # Ensuring that an importance type is selected for plotting
        if self.primary_importance_type is None:
            if importance_type is None:
                raise ValueError("Please enter an importance type.")
            elif importance_type not in self.importance_df.columns:
                raise ValueError("Invalid importance type. Please enter a valid importance type.")
            else:
                self.primary_importance_type = importance_type

        # Sort the dataframe by the primary importance type
        self.importance_df.sort_values(by=self.primary_importance_type, ascending=False, inplace=True)

        # Plotting the feature importance
        fig, ax = plt.subplots()
        ax.barh(
            y=self.importance_df["Bits"][:max_importance],
            width=self.importance_df[self.primary_importance_type][:max_importance],
            color="blue",
        )
        ax.bar_label(ax.containers[0], fmt="%.1f")  # labels on each bar
        ax.margins(x=0.1)  # padding the y-margin so that the numbers fit in properly

        # inverting to avoid ascending order
        plt.gca().invert_yaxis()

        # Hide the left, right and top spines
        ax.spines[["left", "right", "top"]].set_visible(False)

        # Adding the title
        # plt.title(f"Top {max_importance} Bits by XGBoost's {self.primary_importance_type.capitalize()} Feature Importance")

        # Adding labels
        plt.ylabel("")
        plt.xlabel(f"f-score (feature importance by {self.primary_importance_type.capitalize()})")

        # check if save_fig is provided
        if save_fig != "":
            plt.savefig(save_fig, bbox_inches="tight", dpi=300)
            plt.clf()
            plt.close()
            return
        else:
            plt.show()

    def single_feature_importance_plot(self, query_bit: int, return_fig_path: bool = True) -> Optional[str]:
        """
        Saves a feature importance plot for a single bit.

        Parameters
        ----------
        query_bit : int
            Bit for which the feature importance plot is to be generated.
        return_fig_path : bool, optional
            Whether to return the path to the figure or not (default is True).

        Returns
        -------
        str or None
            Path to the figure if `return_fig_path` is True, None otherwise.

        Raises
        ------
        ValueError
            If the user hasn't run the `.make_feature_importance_df()` method first.

        Notes:
        -----
        The plot is saved in the folder with the name:
        Bit_{query_bit}_Feature_Importance_{self.primary_importance_type}_Plot.png
        """
        # validate if self.importance_df exists
        if self.importance_df is None:
            raise ValueError("Please run the make_feature_importance_df() method first.")

        # making the feature importance plot of only the query bit
        fig, ax = plt.subplots(figsize=[8, 1])
        ax.barh(
            y=self.importance_df["Bits"][query_bit],
            width=self.importance_df[self.primary_importance_type][query_bit],
            color="blue",
        )
        ax.bar_label(ax.containers[0], fmt="%.1f")  # labels on each bar
        ax.margins(x=0.1)  # padding the y-margin so that the numbers fit in properly

        # hide the left, right and top spines
        ax.spines[["left", "right", "top"]].set_visible(False)

        # adding the title
        # plt.title(f"Top {max_importance} Bits by XGBoost's {self.primary_importance_type.capitalize()} Feature Importance")

        # adding labels
        plt.ylabel("")
        plt.xlabel(f"f-score (feature importance by {self.primary_importance_type.capitalize()})")

        # save the plot
        feature_importance_plot_path = (
            f"{self.folder_name}/Bit_{query_bit}_Feature_Importance_{self.primary_importance_type}_Plot.png"
        )
        plt.savefig(feature_importance_plot_path, dpi=1000, bbox_inches="tight")

        plt.clf()
        plt.close()

        if return_fig_path:
            return feature_importance_plot_path

    def generate_bit_analysis_report(
        self,
        bit_list: list[int],
        report_title: str = "",
        time_stamp: bool = True,
        file_name: str = "X-FP_Bit_Analysis_Report",
    ) -> None:
        """
        Generates a bit analysis PDF report for the query bits.

        Parameters
        ----------
        bit_list : list[int]
            List of bits for which the bit analysis report is to be generated.
        report_title : str, optional
            Title of the report (default is "").
        time_stamp : bool, optional
            Whether to add a time stamp to the report or not (default is True).
        file_name : str, optional
            Name of the report file (default is "X-FP_Bit_Analysis_Report").

        Raises
        ------
        ValueError
            If the user hasn't run the `.make_feature_importance_df()` method first.
        """

        # folder settings
        self._folder_settings()

        # validate if self.importance_df exists
        if self.importance_df is None:
            raise ValueError("Please run the make_feature_importance_df() method first.")

        # start the report
        bit_analysis_report = BitAnalysis()

        # generate the overall feature importance plot
        self.generate_feaure_importance_plot(
            save_fig=f"{self.folder_name}/Overall_Feature_Importance_{self.primary_importance_type}_Plot.png"
        )
        overall_feature_importance_image = (
            f"{self.folder_name}/Overall_Feature_Importance_{self.primary_importance_type}_Plot.png"
        )

        # print the title page
        bit_analysis_report.intro_page(
            overall_feature_importance_image=overall_feature_importance_image,
            img_text=f"XGBoost's {self.primary_importance_type.capitalize()} Feature Importance Plot for the Top 10 Morgan Fingerprint Bits",
            report_title=report_title,
            time_stamp=time_stamp,
        )

        # looping through the bits
        for query_bit in bit_list:
            # generate the feature importance plot for the query bit, saving plot, and retrieving the path
            single_summary_plot_path = self.single_feature_importance_plot(query_bit=query_bit, return_fig_path=True)

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
                mol_importance_image=single_summary_plot_path,
                mol_image=substructure_images_path,
                in_df=substructure_count_df,
                mol_importance_text=f"XGBoost '{self.primary_importance_type.capitalize()}' feature importance plot for the Bit {query_bit}",
                n_compounds=len(self.input_df),
            )

        # Add end note to the report
        path_to_xgb_feature_importance_note = pkg_resources.files(resources) / "xgb_feature_importance.txt"
        bit_analysis_report.add_end_note(feature_importance_note=path_to_xgb_feature_importance_note)

        # Save the report
        bit_analysis_report.output(f"{self.folder_name}/{file_name}.pdf")
