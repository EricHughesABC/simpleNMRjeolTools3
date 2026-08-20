from qtpy.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QComboBox,
    QPushButton,
    QGridLayout,
    QScrollArea,
    QWidget,
    QFrame,
    QCheckBox,
)
from qtpy.QtCore import Qt


class NMRExperimentDialog(QDialog):
    def __init__(self, experiment_names, chosen_types={}, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Choose and Identify NMR Experiments")
        self.setModal(True)
        self.resize(600, 450)  # Slightly taller for checkboxes

        # Fixed list of experiment types
        self.experiment_types = [
            "SKIP",
            "H1_1D",
            "C13_1D",
            "Pureshift",
            "DEPT135",
            "HSQC",
            "HMBC",
            "COSY",
            "HSQCCLIPCOSY",
            "DDEPTCH3ONLY",
        ]

        # Store dropdowns for later access
        self.experiment_dropdowns = []

        # Set initial values for dropdowns based on chosen_types
        self.chosen_types = chosen_types
        self.experiment_names = experiment_names

        # Create main layout
        main_layout = QVBoxLayout(self)

        # Add title (optional, since it's already in window title)
        title_label = QLabel("Choose and Identify NMR Experiments")
        title_label.setAlignment(Qt.AlignCenter)
        title_label.setStyleSheet("font-size: 14px; font-weight: bold; margin: 10px;")
        main_layout.addWidget(title_label)

        # Create scroll area for experiment rows
        scroll_area = QScrollArea()
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)

        # Create header row
        header_frame = QFrame()
        header_layout = QGridLayout(header_frame)
        header_layout.setContentsMargins(10, 5, 10, 5)

        found_header = QLabel("Found Experiments with Peaks")
        found_header.setStyleSheet("font-weight: bold; font-size: 12px;")
        found_header.setAlignment(Qt.AlignCenter)

        type_header = QLabel("Experiment Type")
        type_header.setStyleSheet("font-weight: bold; font-size: 12px;")
        type_header.setAlignment(Qt.AlignCenter)

        header_layout.addWidget(found_header, 0, 0)
        header_layout.addWidget(type_header, 0, 1)
        header_layout.setColumnStretch(0, 1)
        header_layout.setColumnStretch(1, 1)

        scroll_layout.addWidget(header_frame)

        # Add separator line
        separator = QFrame()
        separator.setFrameShape(QFrame.HLine)
        separator.setFrameShadow(QFrame.Sunken)
        scroll_layout.addWidget(separator)

        # Create experiment rows
        for i, exp_name in enumerate(experiment_names):
            row_frame = QFrame()
            row_layout = QGridLayout(row_frame)
            row_layout.setContentsMargins(10, 5, 10, 5)

            # Left side: Experiment name (static text)
            name_label = QLabel(exp_name)
            name_label.setStyleSheet(
                "padding: 5px; border: 1px solid #ccc; background-color: #f9f9f9;"
            )
            name_label.setAlignment(Qt.AlignCenter)

            # Right side: Experiment type dropdown
            type_dropdown = QComboBox()
            type_dropdown.addItems(self.experiment_types)
            if exp_name in self.chosen_types:
                type_dropdown.setCurrentText(self.chosen_types[exp_name])
            else:
                type_dropdown.setCurrentText("SKIP")  # Default selection
            type_dropdown.setStyleSheet("padding: 5px;")
            self.experiment_dropdowns.append(type_dropdown)

            # Add to grid layout
            row_layout.addWidget(name_label, 0, 0)
            row_layout.addWidget(type_dropdown, 0, 1)
            row_layout.setColumnStretch(0, 1)
            row_layout.setColumnStretch(1, 1)

            scroll_layout.addWidget(row_frame)

        # Add stretch to push everything to top
        scroll_layout.addStretch()

        scroll_area.setWidget(scroll_widget)
        scroll_area.setWidgetResizable(True)
        main_layout.addWidget(scroll_area)

        # Add checkboxes section
        checkbox_frame = QFrame()
        checkbox_frame.setFrameShape(QFrame.Box)
        checkbox_frame.setFrameShadow(QFrame.Raised)
        checkbox_frame.setStyleSheet(
            "QFrame { background-color: #f0f0f0; margin: 5px; padding: 5px; }"
        )

        checkbox_layout = QVBoxLayout(checkbox_frame)

        # Options label
        options_label = QLabel("Processing Options:")
        options_label.setStyleSheet("font-weight: bold; margin-bottom: 5px;")
        checkbox_layout.addWidget(options_label)

        # Checkbox 1: Use Simulated Annealing (True by default)
        self.simulated_annealing_cb = QCheckBox(
            "Use Simulated Annealing to Optimize Assignment"
        )
        self.simulated_annealing_cb.setChecked(True)  # True by default
        self.simulated_annealing_cb.setStyleSheet("margin-left: 10px;")
        checkbox_layout.addWidget(self.simulated_annealing_cb)

        # Checkbox 2: Auto-validate assignments
        self.auto_validate_cb = QCheckBox("Permission to save results to Database")
        self.auto_validate_cb.setChecked(False)  # False by default
        self.auto_validate_cb.setStyleSheet("margin-left: 10px;")
        checkbox_layout.addWidget(self.auto_validate_cb)

        main_layout.addWidget(checkbox_frame)

        # Create button layout
        button_layout = QHBoxLayout()

        # Reset button to set all dropdowns back to SKIP
        reset_btn = QPushButton("Reset All to SKIP")
        reset_btn.clicked.connect(self.reset_all)

        # OK and Cancel buttons
        ok_btn = QPushButton("OK")
        cancel_btn = QPushButton("Cancel")
        ok_btn.clicked.connect(self.accept)
        cancel_btn.clicked.connect(self.reject)

        # Add buttons to layout
        button_layout.addWidget(reset_btn)
        button_layout.addStretch()
        button_layout.addWidget(ok_btn)
        button_layout.addWidget(cancel_btn)

        main_layout.addLayout(button_layout)

    def reset_all(self):
        """Reset all dropdowns to SKIP"""
        for dropdown in self.experiment_dropdowns:
            dropdown.setCurrentText("SKIP")

    def get_experiment_assignments(self):
        """Get the experiment name to type assignments"""
        assignments = []
        for i, dropdown in enumerate(self.experiment_dropdowns):
            assignments.append(
                {
                    "experiment_name": self.experiment_names[i],
                    "experiment_type": dropdown.currentText(),
                    "index": i,
                }
            )
        return assignments

    def get_processing_options(self):
        """Get the processing options from checkboxes"""
        return (
            self.simulated_annealing_cb.isChecked(),
            self.auto_validate_cb.isChecked(),
        )


# def main():
#     """Main function using guidata context"""

#     # Use guidata's Qt application context
#     with qt_app_context():
#         print("Starting NMR Experiment Dialog example...")

#         # Generate random experiment names (simulating found experiments)
#         num_experiments = random.randint(3, 7)
#         experiment_names = []

#         # Create some realistic experiment names (no SKIP in names anymore)
#         experiment_prefixes = ["Experiment", "Sample", "Data", "Spec", "Run"]
#         for i in range(num_experiments):
#             prefix = random.choice(experiment_prefixes)
#             experiment_names.append(f"{prefix} {chr(65 + i)}")

#         print(f"Found {num_experiments} experiments: {experiment_names}")

#         # Create and show the dialog
#         dialog = NMRExperimentDialog(experiment_names)

#         # Execute dialog and get result
#         result = dialog.exec_()

#         if result == QDialog.Accepted:
#             print("\n=== EXPERIMENT ASSIGNMENTS ===")
#             assignments = dialog.get_experiment_assignments()

#             for assignment in assignments:
#                 print(
#                     f"'{assignment['experiment_name']}' -> {assignment['experiment_type']}"
#                 )

#             print("\n=== PROCESSING OPTIONS ===")
#             simulated_annealing, ml_consent = dialog.get_processing_options()
#             print(f"Use Simulated Annealing: {simulated_annealing}")
#             print(f"Save to Database Consent: {ml_consent}")

#         else:
#             print("\nDialog was cancelled by user.")

#         print("\nProgram ending.")


# if __name__ == "__main__":
#     import random
#     main()
