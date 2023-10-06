import importlib.resources as pkg_resources
from datetime import datetime

import pandas as pd
from fpdf import FPDF, Align

from . import resources


class BitAnalysis(FPDF):
    def header(self) -> None:
        """
        Generate header of the report.
        """

        # For the cover page, add the X-FP logo and the CzodrowskiLab logo
        if self.page_no() == 1:
            # X-FP Logo
            self.image(
                pkg_resources.files(resources) / "X-FP_logo.png",
                w=60,
                x=Align.L,
                y=10,
            )

            # CzodrowskiLab Logo
            self.image(
                pkg_resources.files(resources) / "CzodrowskiLab_Logo.png",
                w=30,
                x=Align.R,
                y=5,
            )

        else:
            self.title = "X-FP"

            # Setting header's font:
            self.set_font("helvetica", "B", 12)

            # Calculating width of title
            width = self.get_string_width(self.title) + 6

            # Setting cursor position:
            self.set_x((210 - width) / 2)

            # Setting thickness of the frame (1 mm)
            self.set_line_width(1)
            # Printing title:
            self.cell(
                width,
                9,
                self.title,
                new_x="LMARGIN",
                new_y="NEXT",
                align="C",
                fill=False,
            )

            self.ln(5)

            # printing subtitle
            self.set_font("helvetica", "I", 10)
            self.set_x((210 - width) / 2)
            self.cell(
                width,
                h=5,
                txt="Bit Analysis",
                new_x="LMARGIN",
                new_y="NEXT",
                align="C",
                fill=False,
            )

        # Performing a line break:
        self.ln(10)

    def footer(self) -> None:
        """
        Generate footer of the report.
        """
        # Setting position at 1.5 cm from bottom:
        self.set_y(-15)
        # Setting font: helvetica italic 8
        self.set_font("helvetica", "I", 8)
        # Setting text color to gray:
        self.set_text_color(128)
        # Adding the user provided report title on the left side of the footer
        self.cell(0, 10, self.report_title, align="L")

        # Setting position at 1.5 cm from bottom:
        self.set_y(-15)
        # Setting font: helvetica italic 8
        self.set_font("helvetica", "I", 8)
        # Setting text color to gray:
        self.set_text_color(128)
        # Printing page number
        self.cell(0, 10, f"Page {self.page_no()} of {{nb}}", align="C")

        # time stamp
        if self.time_stamp:
            # Setting position at 1.5 cm from bottom:
            self.set_y(-15)
            # Setting font: helvetica italic 8
            self.set_font("helvetica", "I", 8)
            # Setting text color to gray:
            self.set_text_color(128)
            # Printing time stamp
            self.cell(0, 10, f"Report generated on {self.time_stamp}", align="R")

    def intro_page(
        self, overall_feature_importance_image: str, img_text: str, report_title: str = "", time_stamp: bool = True
    ) -> None:
        """
        Generate an introduction page of the report where the overall feature importance plot is displayed.

        Parameters
        ----------
        overall_feature_importance_image : str
            Image of the overall feature importance plot.
        img_text : str
            Text to be displayed above the overall feature importance image.
        report_title : str, optional
            Title of the report provided by the user (default is "").
        time_stamp : bool, optional
            Add time stamp to the report footer or not (default is True). The time is local and is formatted as
            '%Y-%m-%d %H:%M:%S'.
        """

        # Truncate the title length to 50 characters
        if len(report_title) > 50:
            report_title = f"{report_title[:50]}..."

        # Creating a class variable to store the report title
        self.report_title = report_title

        # Creating a class variable to store the time stamp
        if time_stamp:
            self.time_stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        else:
            self.time_stamp = None

        # Add page
        self.add_page()
        self.set_y(70)  # in case report title is missing

        # writing the title of the report if provided
        if report_title != "":
            self.set_font("helvetica", "U", 20)
            self.set_y(40)
            self.multi_cell(w=0, txt=self.report_title, align="C")
            self.ln(20)

        self.image(
            overall_feature_importance_image,
            w=150,
            x=Align.C,
        )
        self.ln(5)

        # adding the title of the plot
        self.set_font("helvetica", "", 10)
        self.multi_cell(0, 5, img_text, align="C")
        self.ln(20)

    def chapter_title(self, in_bit: int) -> None:
        """
        Generate chapter title of the report.

        Parameters
        ----------
        in_bit : int
            Bit of choice for which the analysis is to be performed
        """

        # Setting font: helvetica 12
        self.set_font("helvetica", "", 12)
        # Setting background color
        self.set_fill_color(200, 220, 255)
        # Printing chapter name:
        self.cell(
            0,
            6,
            # f"Chapter {num} : {label}",
            f"Bit {in_bit}",
            new_x="LMARGIN",
            new_y="NEXT",
            align="L",
            fill=True,
        )
        # Performing a line break:
        self.ln(4)

    def print_chapter(
        self,
        in_bit: int,
        mol_importance_image: str,
        mol_image: str,
        in_df: pd.DataFrame,
        mol_importance_text: str = None,
        n_compounds: int = None,
    ) -> None:
        """
        Generate chapter of the report.

        Parameters
        ----------
        in_bit : int
            Bit of choice for which the analysis is to be performed.
        mol_importance_image : str
            Image of the feature importance plot for the query bit.
        mol_image : str
            Image of the substructure plot for the query bit.
        in_df : pd.DataFrame
            Pandas DataFrame containing the substructure frequency table for the query bit.
        mol_importance_text : str, optional
            Text to be displayed above the feature importance image (default is None).
        n_compounds : int, optional
            Total number of compounds in the input dataframe (default is None).
        """

        # Setting font: helvetica 8
        self.set_font("helvetica", "", 8)

        self.add_page()
        self.chapter_title(in_bit)

        # Feature importance image, e.g., SHAP Plot

        if mol_importance_text is None:
            mol_importance_text = f"SHAP Summary Plot for the Bit {in_bit}:"

        self.multi_cell(0, 5, mol_importance_text)
        self.ln()
        self.image(mol_importance_image, w=100, x=20)
        self.ln()

        # Substructure images
        sub_struct_heading_text = f"Substructures present:"
        self.multi_cell(0, 5, sub_struct_heading_text)
        self.ln()
        self.image(mol_image, w=150, x=20)
        self.ln()

        # Substructure table

        if n_compounds is None:
            sub_struct_table_text = "Substructures frequency table:"
        else:
            sub_struct_table_text = f"Substructures frequency table (total compounds: {n_compounds}):"
        self.multi_cell(0, 5, sub_struct_table_text)
        self.ln()

        in_df = in_df.copy()
        in_df = in_df.applymap(str)  # Convert all data inside dataframe into string type

        columns = [list(in_df)]  # Get list of dataframe columns
        rows = in_df.values.tolist()  # Get list of dataframe rows
        data = columns + rows  # Combine columns and rows in one list

        self.set_font("helvetica", size=8)
        with self.table(
            borders_layout="HORIZONTAL_LINES",
            # cell_fill_color=200,  # grey
            cell_fill_mode="ROWS",
            line_height=self.font_size * 2.5,
            text_align="CENTER",
            width=180,
            col_widths=(20, 100, 30, 30),
        ) as table:
            for data_row in data:
                row = table.row()
                for datum in data_row:
                    row.cell(datum)
        self.ln()

    def add_end_note(self, feature_importance_note: str = None) -> None:
        """
        Generate end note of the report.

        Parameters
        ----------
        feature_importance_note : str, optional
            Path to the text file containing the information about feature importance (default is None).
        """

        # Add page
        self.add_page()

        # setting font of title
        self.set_font("helvetica", "B", 12)

        self.cell(
            0,
            6,
            f"Notes about substructures rendering: ",
            new_x="LMARGIN",
            new_y="NEXT",
            align="L",
            # fill=True,
        )

        # setting font of text
        self.set_font("helvetica", "", 8)

        # Reading text file of substructure rendering:
        substructure_rendering_file = pkg_resources.files(resources) / "rendering_info.txt"
        with substructure_rendering_file.open() as fh:
            substructure_rendering_text = fh.read()  # .decode("latin-1")

        self.multi_cell(0, 5, substructure_rendering_text)
        self.ln()

        # Adding note about feature importance
        if feature_importance_note:
            self.set_font("helvetica", "B", 12)

            self.cell(
                0,
                6,
                f"Notes about feature importance: ",
                new_x="LMARGIN",
                new_y="NEXT",
                align="L",
                # fill=True,
            )

            # setting font of text
            self.set_font("helvetica", "", 8)

            # Reading text file of feature importance:
            if feature_importance_note == "shap_tree_explainer":
                feature_importance_file = pkg_resources.files(resources) / "shap_tree_explainer_info.txt"

            else:
                feature_importance_file = feature_importance_note

            with feature_importance_file.open() as fh:
                feature_importance_text = fh.read()

            self.multi_cell(0, 5, feature_importance_text)
            self.ln()
