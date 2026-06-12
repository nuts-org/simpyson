# src/simpyson/gui.py
from __future__ import annotations

import re
import sys
import tempfile
from pathlib import Path

import plotly.graph_objects as go
from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt, QUrl
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtWidgets import (
    QAction,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QSplitter,
    QVBoxLayout,
    QWidget,
)

from simpyson.io import read_simp
from simpyson.utils import add_spectra, get_larmor_freq


class SimpysonGUI(QMainWindow):
    def __init__(self, files_to_open=None):
        super().__init__()
        self.setWindowTitle("Simpyson GUI")
        self.setGeometry(100, 100, 1200, 800)

        self.files_data = {}
        self.current_file = None
        self.data = None
        self.filename = None
        self.b0_edit = None
        self.nucleus_edit = None
        self.view_combo = None

        self.init_ui()

        if files_to_open:
            self.open_files(files_to_open)

    def init_ui(self):
        self.create_menu()
        self.create_main_layout()

    def file_list_key_press(self, event):
        if event.key() == Qt.Key_Delete:
            selected_items = self.file_list.selectedItems()
            if selected_items:
                confirmation = QMessageBox.question(
                    self, 'Remove Files',
                    f'Remove {len(selected_items)} selected file(s)?',
                    QMessageBox.Yes | QMessageBox.No
                )

                if confirmation == QMessageBox.Yes:
                    for item in selected_items:
                        file_name = item.text()

                        if file_name in self.files_data:
                            del self.files_data[file_name]

                        self.file_list.takeItem(self.file_list.row(item))

                    if not self.files_data:
                        self.current_file = None
                        self.data = None
                        self.filename = None
                    elif self.current_file not in self.files_data:
                        self.current_file = next(iter(self.files_data))
                        self.data = self.files_data[self.current_file]['data']
                        self.filename = self.files_data[self.current_file]['path']

                    self.plot_data()
        else:
            QListWidget.keyPressEvent(self.file_list, event)

    def create_main_layout(self):
        # Create central widget with horizontal layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QHBoxLayout(central_widget)

        # Create splitter for resizable panels
        splitter = QSplitter(Qt.Horizontal)

        # LEFT PANEL: file list + selection panel
        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)
        left_layout.setContentsMargins(0, 0, 0, 0)
        left_layout.setSpacing(6)

        # File list
        self.file_list = QListWidget()
        self.file_list.setSelectionMode(QListWidget.ExtendedSelection)
        self.file_list.itemSelectionChanged.connect(self.on_selection_changed)
        self.file_list.keyPressEvent = self.file_list_key_press
        left_layout.addWidget(self.file_list)

        # Selection settings panel
        controls = QGroupBox('Selection settings')
        form = QFormLayout(controls)
        form.setContentsMargins(8, 8, 8, 8)

        self.b0_edit = QLineEdit()
        self.b0_edit.setPlaceholderText('e.g., 400MHz or 9.4T')
        form.addRow(QLabel('B0:'), self.b0_edit)

        self.nucleus_edit = QLineEdit()
        self.nucleus_edit.setPlaceholderText('e.g., 1H or 13C')
        form.addRow(QLabel('Nucleus:'), self.nucleus_edit)

        # View toggle (Hz/ppm)
        self.view_combo = QComboBox()
        self.view_combo.addItems(['Hz', 'ppm'])
        form.addRow(QLabel('View:'), self.view_combo)

        # Action buttons
        btn_row = QHBoxLayout()
        apply_btn = QPushButton('Apply to selection')
        apply_btn.clicked.connect(self.apply_selection_settings)
        set_view_btn = QPushButton('Set view for selection')
        set_view_btn.clicked.connect(self.change_view_from_combo)
        btn_row.addWidget(apply_btn)
        btn_row.addWidget(set_view_btn)
        form.addRow(btn_row)

        left_layout.addWidget(controls)

        splitter.addWidget(left_panel)

        # RIGHT PANEL: plot area
        plot_widget = QWidget()
        plot_layout = QVBoxLayout(plot_widget)
        self.browser = QWebEngineView()
        plot_layout.addWidget(self.browser)
        plot_widget.setLayout(plot_layout)

        splitter.addWidget(plot_widget)
        splitter.setSizes([260, 940])  # Left panel width, rest for plot

        layout.addWidget(splitter)

    def create_menu(self):
        menubar = self.menuBar()
        menubar.setNativeMenuBar(False)

        # File Menu
        file_menu = menubar.addMenu('File')

        open_action = QAction('Open', self)
        open_action.triggered.connect(self.open_file)
        open_action.setShortcut('Ctrl+o')
        file_menu.addAction(open_action)

        save_action = QAction('Save', self)
        save_action.triggered.connect(self.save_file)
        save_action.setShortcut('Ctrl+s')
        file_menu.addAction(save_action)

        file_menu.addSeparator()

        exit_action = QAction('Exit', self)
        exit_action.triggered.connect(QtWidgets.qApp.quit)
        file_menu.addAction(exit_action)

        # Process Menu
        process_menu = menubar.addMenu('Process')

        setup_conversions = QAction('Setup Conversions', self)
        setup_conversions.triggered.connect(self.setup_conversions)
        setup_conversions.setShortcut('Ctrl+g')
        process_menu.addAction(setup_conversions)

        # Combine spectra menu
        combine_spectra = QAction('Combine Spectra', self)
        combine_spectra.triggered.connect(self.combine_selected_spectra)
        process_menu.addAction(combine_spectra)

        # View Menu
        view_menu = menubar.addMenu('View')

        view_hz = QAction('Hz', self)
        view_hz.triggered.connect(lambda: self.change_view('hz'))
        view_hz.setShortcut('Ctrl+1')
        view_menu.addAction(view_hz)

        view_ppm = QAction('ppm', self)
        view_ppm.triggered.connect(lambda: self.change_view('ppm'))
        view_ppm.setShortcut('Ctrl+2')
        view_menu.addAction(view_ppm)

        view_fid = QAction('FID', self)
        view_fid.triggered.connect(lambda: self.change_view('fid'))
        view_fid.setShortcut('Ctrl+3')
        view_menu.addAction(view_fid)

    def save_file(self):
        if not self.current_file:
            QMessageBox.warning(self, 'Save File', 'No file selected!')
            return

        if not self.data:
            QMessageBox.warning(self, 'Save File', 'No data to save!')
            return

        options = QFileDialog.Options()
        save_filename, _ = QFileDialog.getSaveFileName(
            self,
            'Save File',
            '',
            'All Supported Files (*.csv *.fid *.spe *.xreim *.csdf);;CSV Files (*.csv);;'
            'FID Files (*.fid);;SPE Files (*.spe);;XREIM Files (*.xreim);;CSDF Files (*.csdf)',
            options=options
        )

        if not save_filename:
            return

        try:
            format = save_filename.lower().split('.')[-1]

            if format not in ['spe', 'fid', 'csv', 'xreim', 'csdf']:
                QMessageBox.warning(self, 'Save File', 'Unsupported file format!')
                return

            self.data.write(save_filename, format=format)

            QMessageBox.information(self, 'Save File', 'File saved successfully!')
        except Exception as e:
            QMessageBox.warning(self, 'Save File', f'Error saving file: {e!s}')

    def open_files(self, filenames):
        for filename in filenames:
            if filename:
                file_format = Path(filename).suffix.lstrip('.').lower()
                base_name = Path(filename).name

                if file_format in ('spe', 'csdf'):
                    view = 'hz'
                elif file_format == 'fid':
                    view = 'fid'
                elif file_format == 'xreim':
                    view = 'xreim'
                else:
                    QMessageBox.warning(
                        self, 'Open File',
                        f'Unsupported file format: {base_name}. '
                        'Supported: .spe, .fid, .xreim, .csdf',
                    )
                    continue

                try:
                    data = read_simp(filename, format=file_format)
                except (OSError, ValueError) as e:
                    QMessageBox.warning(self, 'Open File', f'Error reading {base_name}: {e!s}')
                    continue

                self.files_data[base_name] = {
                    'data': data,
                    'path': filename,
                    'view': view,
                }

                # Add to list widget if not already there
                if self.file_list.findItems(base_name, Qt.MatchExactly) == []:
                    self.file_list.addItem(base_name)

                # If this is the first file, display it
                if not self.current_file:
                    self.current_file = base_name
                    self.data = self.files_data[base_name]['data']
                    self.filename = filename

        if self.current_file:
            self.plot_data()

    def open_file(self):
        options = QFileDialog.Options()
        filenames, _ = QFileDialog.getOpenFileNames(
            self, 'Open File', '', 'SIMPSON Files (*.spe *.fid *.xreim *.csdf)', options=options
        )
        if filenames:
            self.open_files(filenames)

    def on_selection_changed(self):
        selected_items = self.file_list.selectedItems()
        if selected_items:
            self.current_file = selected_items[0].text()
            self.data = self.files_data[self.current_file]['data']
            self.filename = self.files_data[self.current_file]['path']
            self.plot_data(selected_items)

    def plot_data(self, selected_items=None):
        if not selected_items:
            selected_items = [self.file_list.findItems(self.current_file, Qt.MatchExactly)[0]] if self.current_file else []

        if selected_items:
            fig = go.Figure()

            # Determine x-axis type from first selected item
            first_file = self.files_data[selected_items[0].text()]
            first_data = first_file['data']

            if 'view' not in first_file:
                QMessageBox.warning(self, 'Plot Data', 'Data does not contain valid axis information.')
                return

            # Get data based on view
            match first_file['view']:
                case 'ppm':
                    if not first_data.ppm:
                        QMessageBox.warning(self, 'Plot Data', 'No ppm data available. Set B0 and nucleus first.')
                        return
                    x_axis = 'ppm'
                    data_source = 'ppm'
                    xlabel = 'Chemical Shift (ppm)'
                case 'hz':
                    if not first_data.spe:
                        QMessageBox.warning(self, 'Plot Data', 'No spectrum data available.')
                        return
                    x_axis = 'hz'
                    data_source = 'spe'
                    xlabel = 'Frequency (Hz)'
                case 'fid':
                    if not first_data.fid:
                        QMessageBox.warning(self, 'Plot Data', 'No FID data available.')
                        return
                    x_axis = 'time'
                    data_source = 'fid'
                    xlabel = 'Time (ms)'
                case 'xreim':
                    if not first_data.xreim:
                        QMessageBox.warning(self, 'Plot Data', 'No xreim data available.')
                        return
                    x_axis = 'time'
                    data_source = 'xreim'
                    xlabel = 'Time'

            # Plot each selected spectrum
            for item in selected_items:
                data = self.files_data[item.text()]['data']

                # Access data via property (fid/spe/ppm)
                data_dict = getattr(data, data_source)
                if data_dict and x_axis in data_dict:
                    fig.add_trace(go.Scatter(
                        x=data_dict[x_axis],
                        y=data_dict['real'],
                        name=item.text(),
                        mode='lines'
                    ))
                else:
                    QMessageBox.warning(self, 'Plot Data', f'Missing data for {item.text()}.')
                    continue

            # Update layout with consistent axis settings
            fig.update_layout(
                xaxis_title=xlabel,
                yaxis_title='Intensity',
                showlegend=True
            )

            # Remove grids
            fig.update_xaxes(showgrid=False, showline=True, linewidth=2, linecolor='black', ticks='outside', tickwidth=2, tickcolor='black')
            fig.update_yaxes(showgrid=False, showline=True, linewidth=2, linecolor='black', ticks='outside', tickwidth=2, tickcolor='black')

            fig.update_layout(
                plot_bgcolor='rgba(0,0,0,0)',
                paper_bgcolor='rgba(0,0,0,0)'
            )

            # Always invert x-axis for ppm and hz
            if first_file['view'] in ['ppm', 'hz']:
                fig.update_xaxes(autorange="reversed")

            # Update plot display
            html_content = fig.to_html(include_plotlyjs=True, full_html=True)

            html_content = re.sub(r':focus-visible\s*\{[^}]*\}', '', html_content)

            # Save to temp file and load
            temp_path = Path(tempfile.gettempdir()) / 'simpyson_plot.html'
            with temp_path.open('w', encoding='utf-8') as f:
                f.write(html_content)

            self.browser.load(QUrl.fromLocalFile(temp_path))
        else:
            self.browser.setHtml('<html><body><h3 style="text-align:center;margin-top:40px;color:#888;">No data to display</h3></body></html>')
            QMessageBox.warning(self, 'Plot Data', 'No data to plot!')

    def apply_selection_settings(self):
        if not (selected_items := self.get_selection()):
            return

        b0 = self.b0_edit.text().strip()
        nucleus = self.nucleus_edit.text().strip()

        if not (b0 and nucleus):
            QMessageBox.warning(self, 'Selection settings', 'Please fill both B0 and Nucleus.')
            return

        try:
            # Validate inputs
            _ = get_larmor_freq(b0, nucleus)
        except ValueError as e:
            QMessageBox.warning(self, 'Selection settings', f'Invalid input: {e}')
            return

        # Apply to all selected files
        for item in selected_items:
            file_name = item.text()
            data = self.files_data[file_name]['data']
            data.b0 = b0
            data.nucleus = nucleus

        # If current view is ppm (selected in combo), update view too
        if self.view_combo.currentText().lower() == 'ppm':
            self.change_view_from_combo()
        else:
            self.plot_data(selected_items)

    # Change view based on combo for selected items
    def change_view_from_combo(self):
        if not (selected_items := self.get_selection()):
            return

        desired = self.view_combo.currentText().lower()  # 'hz' or 'ppm'
        if desired == 'ppm' and not all(self.has_setup(item.text()) for item in selected_items):
            QMessageBox.warning(self, 'Selection settings', 'B0 and Nucleus must be set for all selected items to use ppm.')
            return

        for item in selected_items:
            file_name = item.text()
            self.files_data[file_name]['view'] = desired

        self.plot_data(selected_items)

    def change_view(self, view):
        if not (selected_items := self.get_selection()):
            return

        if view == 'ppm' and not all(self.has_setup(item.text()) for item in selected_items):
            QMessageBox.warning(self, 'Not Setup', 'B0 and nucleus must be set for PPM view')
            return

        for item in selected_items:
            file_name = item.text()
            self.files_data[file_name]['view'] = view

        self.plot_data(selected_items)

    def setup_conversions(self):
        if not (selected_items := self.get_selection()):
            return

        dialog = QDialog(self)
        dialog.setWindowTitle('Setup Conversions')
        form = QFormLayout(dialog)

        b0_edit = QLineEdit(dialog)
        form.addRow('B0 (e.g, 400MHz or 9.4T:', b0_edit)

        nucleus_edit = QLineEdit(dialog)
        form.addRow('Nucleus (e.g, 1H or 13C:', nucleus_edit)

        buttons = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel,
            parent=dialog
        )
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        form.addRow(buttons)

        if dialog.exec_() == QDialog.Accepted:
            b0 = b0_edit.text()
            nucleus = nucleus_edit.text()

            if not (b0 and nucleus):
                return

            try:
                _ = get_larmor_freq(b0, nucleus)
            except ValueError as e:
                QMessageBox.warning(self, 'Setup Conversions', f'Invalid Input: {e!s}')
                return

            for item in selected_items:
                file_name = item.text()
                file_data = self.files_data[file_name]['data']

                file_data.b0 = b0
                file_data.nucleus = nucleus

    # Combine spectra
    def combine_selected_spectra(self):
        """Combine multiple selected spectra into a single spectrum."""
        if not (selected_items := self.get_selection()):
            return

        if len(selected_items) < 2:
            QMessageBox.warning(self, 'Combine Spectra', 'Please select at least two spectra to combine')
            return

        # Check that all selected items have spectrum data
        for item in selected_items:
            file_data = self.files_data[item.text()]['data']
            if not file_data.spe:
                QMessageBox.warning(self, 'Combine Spectra',
                               f'{item.text()} is not in spectrum format. Convert all files to spectra first.')
                return

        # Collect all Simpy objects
        spectra_list = [self.files_data[item.text()]['data'] for item in selected_items]

        try:
            # Combine spectra
            combined = add_spectra(spectra_list)

            # Create a new name for the combined spectrum
            base_names = [item.text().split('.')[0] for item in selected_items]
            new_name = f"Combined_{'_'.join(base_names[:2])}"
            if len(base_names) > 2:
                new_name += f"_plus{len(base_names)-2}"
            new_name += '.spe'

            # Add to file list
            self.files_data[new_name] = {
                'data': combined,
                'path': None,
                'view': 'hz'
            }

            self.file_list.addItem(new_name)
            new_item = self.file_list.findItems(new_name, Qt.MatchExactly)[0]
            self.file_list.setCurrentItem(new_item)

            QMessageBox.information(self, 'Combine Spectra',
                               f'Successfully combined {len(spectra_list)} spectra')

        except Exception as e:
            QMessageBox.critical(self, 'Error', f'Failed to combine spectra: {e!s}')

    def get_selection(self):
        selected_items = self.file_list.selectedItems()

        if not selected_items:
            QMessageBox.warning(self, 'No Selection', 'Please select at least one item')
            return None

        return selected_items

    def has_setup(self, file_name):
        data = self.files_data[file_name]['data']
        return not (data.b0 is None or data.nucleus is None)

def main(files=None):
    app = QtWidgets.QApplication(sys.argv)
    gui = SimpysonGUI(files_to_open=files)
    gui.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
