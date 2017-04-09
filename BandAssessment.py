import sys
import os
import dicom
import numpy as np
from PIL import Image
from PIL import ImageDraw
from PIL import ImageFilter
import matplotlib.pyplot as plt
import logging
import sqlite3

# define the logging config, output in file
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    filename='running.log',
                    filemode='w')
# define a stream that will show log level > Warning on screen also
console = logging.StreamHandler()
console.setLevel(logging.WARNING)
formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
console.setFormatter(formatter)
logging.getLogger('').addHandler(console)


##############################################################


class SQL3Handler:
    """
    Create a Sqlite3 handler to store the data.
    for the integration result, it will be converted to string and then store into database.
    So when extracting data, the string should be converted back to list or numpy array before doing calculation
    """
    Database_Name = "BandAssessment.sqlite3.db"

    def __init__(self, name, kvp, current, kernel, total_col, slice_thick, instance, integration):
        self.Dicom_Station_Name = name
        self.Dicom_KVP = kvp
        self.Dicom_Current = current
        self.Dicom_Kernel = kernel
        self.Dicom_Total_Collimation = total_col
        self.Dicom_Slice_Thickness = slice_thick
        self.Dicom_Instance = instance
        self.Integration_Result = integration
        logging.debug(r"Run into SQL3Handler")
        try:
            con = sqlite3.connect(self.Database_Name)
        except sqlite3.Error as e:
            logging.debug(str(e))
            return
        logging.debug(r"Database connected")
        sql_cursor = con.cursor()
        sql_string = '''create table BandAssessment(
                           uid integer primary key autoincrement,
                           serial_number integer,
                           tube_voltage real,
                           tube_current integer,
                           kernel text,
                           total_collimation real,
                           slice_thickness real,
                           instance integer,
                           integration_result text);'''
        try:
            sql_cursor.execute(sql_string)
        except sqlite3.Error as e:
            logging.debug(str(e))
            return
        logging.debug(r"create table done.")
        con.close()
    # End of __ini__
    ############################################

    def insert_data(self):
        try:
            con = sqlite3.connect(self.Database_Name)
        except sqlite3.Error as e:
            logging.debug(str(e))
            return
        # convert numpy into string to store in sqlite3
        integration_result = []
        for x in self.Integration_Result:
            integration_result.append(str(x))
        int_result_string = ';'.join(integration_result)
        # set up for store in sql
        sql_cursor = con.cursor()
        sql_string = r"insert into BandAssessment values (?,?,?,?,?,?,?,?,?);"
        try:
            sql_cursor.execute(sql_string,
                               (None, self.Dicom_Station_Name, self.Dicom_KVP, self.Dicom_Current, self.Dicom_Kernel,
                                self.Dicom_Total_Collimation, self.Dicom_Slice_Thickness, self.Dicom_Instance,
                                int_result_string))
        except sqlite3.Error as e:
            logging.error(str(e))
            con.close()
            return
        con.commit()
        logging.info(r"Insert record done.")
        con.close()
    # end of insertExample()
    ############################################

    def read_data(self):
        try:
            con = sqlite3.connect(self.Database_Name)
        except sqlite3.Error as e:
            logging.debug(str(e))
            return
        sql_cursor = con.cursor()
        sql_string = r"select integration_result from BandAssessment"
        sql_cursor.execute(sql_string)
        data = sql_cursor.fetchone()[0]
        str_result = data.split(';')
        float_result = []
        for x in str_result:
            float_result.append(float(x))
        np_result = np.array(float_result)
        print(type(np_result))
        print(np_result)
        con.close()
# End of class SQL3Handler:
##############################################################


class Dicom:
    """
    This is the main class to deal with single dicom image.
    And it will calculate the circular integration of the image.
    Finally the low pass the integration to get a feeling that how the image quality is.
    Function description:
    __imi__: to initialize Dicom file and call other functions to calculate circular integration
    calc_circle: to find the phantom center and find the radius. Normally it has 2 kinds of phantom,
                 20cm and 30cm (diameter)
    bresenham: draw a circle on the image with a radius. return the sum of the points that on the edge of the circle
    integration: this function do 2 parts:
        part 1: call bresenham with different radius (from 1 to phantom radius). And then calculate the
                integration. The integration must be weighted average number because when radius is bigger,
                the number of points is bigger
        part 2: low pass the integration result to more visible.
    """

    def __init__(self, filename, center=0, width=100):
        # set up some basic date
        self.isShowImgReady = False
        self.Center_Col = 256
        self.Center_Row = 256
        self.Radius = 200
        self.Dicom_File_Name = filename
        # open the dicom file
        logging.info(r"Opening file:" + filename)
        try:
            self.Dicom_File = dicom.read_file(filename)
        except Exception as e:
            logging.error(str(e))
            return
        # if file is opened, continue to extract data from dicom file
        try:
            self.Dicom_Station_Name = self.Dicom_File[0x0018, 0x1000].value
            logging.debug(str(self.Dicom_Station_Name) + r"<==System Serial No")
            self.StudyDescription = self.Dicom_File[0x0008, 0x1030].value
            if self.StudyDescription != r"Band Assessment":
                logging.error(self.Dicom_File_Name + " is not Band Assessment")
                return
            logging.debug(self.Dicom_File_Name + ":" + self.StudyDescription + r"<==Study Description is:")
            self.Slop = self.Dicom_File[0x0028, 0x1053].value
            self.Intercept = self.Dicom_File[0x0028, 0x1052].value
            self.Dicom_Image_Data = np.array(self.Dicom_File.pixel_array)
            self.Dicom_Rows = self.Dicom_File[0x0028, 0x0010].value
            self.Dicom_Cols = self.Dicom_File[0x0028, 0x0011].value
            self.Dicom_Pix_Space = self.Dicom_File[0x0028, 0x0030].value
            self.Dicom_KVP = self.Dicom_File[0x0018, 0x0060].value
            self.Dicom_Current = self.Dicom_File[0x0018, 0x1151].value
            self.Dicom_Kernel = self.Dicom_File[0x0018, 0x1210].value
            self.Dicom_Series = self.Dicom_File[0x0020, 0x0011].value
            self.Dicom_Total_Collimation = self.Dicom_File[0x0018, 0x9307].value
            self.Dicom_Slice_Thickness = self.Dicom_File[0x0018, 0x0050].value
            self.Dicom_Instance = self.Dicom_File[0x0020, 0x0013].value
            # self.Dicom_Image_No = int(self.Dicom_Total_Collimation / self.Dicom_Slice_Thickness)
            logging.debug(r"Image Mode: " +
                          str(self.Dicom_KVP) + r"KV_" +
                          str(self.Dicom_Current) + r"mA_" +
                          str(self.Dicom_Kernel) + r"_" +
                          str(self.Dicom_Total_Collimation) + r"I" +
                          str(self.Dicom_Slice_Thickness) + r"." +
                          str(self.Dicom_Instance))
            self.Dicom_Pix_Space = self.Dicom_File[0x0028, 0x0030].value
        except Exception as e:
            logging.error(e)
            return
        # Do the initial calculation
        self.Dicom_HU_Image = self.Dicom_Image_Data * self.Slop + self.Intercept  # Convert to HU unit
        self.Window_Upper = center + width / 2
        self.Window_Lower = center - width / 2
        # make filter according to center and width
        upper_filter = self.Dicom_HU_Image > self.Window_Upper
        self.Dicom_HU_Image[upper_filter] = self.Window_Upper  # set upper value
        lower_filter = self.Dicom_HU_Image < self.Window_Lower
        self.Dicom_HU_Image[lower_filter] = self.Window_Lower  # set lower value
        # rescale the data to 0~255
        min_hu_image = self.Dicom_HU_Image.min()
        self.Image_rescale = self.Dicom_HU_Image + (0 - min_hu_image)  # get rid of minus number
        max_image_rescale = self.Image_rescale.max()
        self.Image_rescale = self.Image_rescale * 255 / max_image_rescale  # rescale the image to fit 0~255
        # make a copy of HU image for calculate center of the circle
        self.Dicom_HU_Copy = self.Dicom_HU_Image.copy()
        # try to calculate radius and center col / row
        self.calc_circle()
        logging.debug(r"Center of circle has been found.")
        # define circular integration result
        self.Integration_Result = np.zeros(self.Radius)
        self.Median_Filter_Result = np.zeros(self.Radius)
        # main calculation
        self.integration()
        logging.debug(r"Circular integration done.")
        self.isShowImgReady = True
    ###############################################

    def calc_circle(self):
        # remove the pixel if the pixel is low
        remove_rate = 1.2  # define rate
        remove_low_value = self.Dicom_HU_Copy < self.Window_Upper / remove_rate  # set up the filter
        self.Dicom_HU_Copy[remove_low_value] = 0  # remove the low values
        # convert the np array to image, then use PIL to find image edge
        im = Image.fromarray(self.Dicom_HU_Copy).convert("L").filter(ImageFilter.FIND_EDGES)
        # convert image back to np array
        filtered_image = np.array(im)
        # use filtered image to re-calculate center and radius
        # the method is simple
        # from up/down/left/right side to go into center
        # the 1st number is > mean value, it's the edge
        # calculate the distance from th edge to center
        # start to calculate center col
        abnormal = False  # set up the flag to store if calculation is abnormal
        left_distance = 0
        right_distance = 0
        up_distance = 0
        low_distance = 0
        for left_distance in range(1, im.size[1]):
            if filtered_image[self.Center_Row, left_distance] != 0:
                break
        for right_distance in range(1, im.size[1]):
            if filtered_image[self.Center_Row, im.size[1] - right_distance] != 0:
                break
        self.Center_Col += (left_distance - right_distance) // 2
        logging.debug(r"Center Col calculated as: " + str(self.Center_Col))
        # if the calculated center col deviated too much
        deviation = 20
        if (self.Center_Col > self.Dicom_Cols // 2 + deviation) or (self.Center_Col < self.Dicom_Cols // 2 - deviation):
            logging.warning(r"It seems abnormal when calculate Center Col, use image center now!")
            self.Center_Col = self.Dicom_Cols // 2
            abnormal = True
        # start to calculate center row
        for up_distance in range(1, im.size[0]):
            if filtered_image[up_distance, self.Center_Col] != 0:
                break
        for low_distance in range(1, im.size[1]):
            if filtered_image[im.size[1] - low_distance, self.Center_Col] != 0:
                break
        self.Center_Row += (up_distance - low_distance) // 2
        logging.debug(r"Center Row calculated as: " + str(self.Center_Row))
        # if the calculated center row deviated too much
        if (self.Center_Row > self.Dicom_Rows // 2 + deviation) or (self.Center_Row < self.Dicom_Rows // 2 - deviation):
            logging.warning(r"It seems abnormal when calculate Center row, use image center now!")
            self.Center_Row = self.Dicom_Rows // 2
            abnormal = True
        # set different radius according to normal/abnormal situation
        if abnormal is False:
            self.Radius = (im.size[0] - left_distance - right_distance) // 2
            diameter_in_cm = self.Radius * self.Dicom_Pix_Space[0] * 2
            logging.debug(str(self.Radius) + r"pix (radius), " + str(diameter_in_cm) +
                          r"cm(diameter)<==Calculated phantom diameter")
            # standardize the radius
            if diameter_in_cm < 250:
                self.Radius = 233
                logging.debug(str(self.Radius) + r"pix" + r", which is: " +
                              str(self.Radius * self.Dicom_Pix_Space[0] * 2) + r"cm <=========Radius Readjusted")
            else:
                self.Radius = 220
                logging.debug(str(self.Radius) + r"pix" + r", which is: " +
                              str(self.Radius * self.Dicom_Pix_Space[0] * 2) + r"cm <=========Radius Readjusted")
        else:
            logging.warning(r"Calculated center is abnormal, use 50 as radius!")
            self.Radius = 50
    #######################################

    def bresenham(self, radius):
        x = 0
        y = radius
        d = 3 - 2 * radius
        while x < y:
            self.Integration_Result[radius] = self.Integration_Result[radius] + self.Dicom_HU_Image[
                self.Center_Row - y, self.Center_Col + x]
            self.Integration_Result[radius] = self.Integration_Result[radius] + self.Dicom_HU_Image[
                self.Center_Row + y, self.Center_Col + x]
            self.Integration_Result[radius] = self.Integration_Result[radius] + self.Dicom_HU_Image[
                self.Center_Row - y, self.Center_Col - x]
            self.Integration_Result[radius] = self.Integration_Result[radius] + self.Dicom_HU_Image[
                self.Center_Row + y, self.Center_Col - x]
            self.Integration_Result[radius] = self.Integration_Result[radius] + self.Dicom_HU_Image[
                self.Center_Row - x, self.Center_Col + y]
            self.Integration_Result[radius] = self.Integration_Result[radius] + self.Dicom_HU_Image[
                self.Center_Row - x, self.Center_Col - y]
            self.Integration_Result[radius] = self.Integration_Result[radius] + self.Dicom_HU_Image[
                self.Center_Row + x, self.Center_Col + y]
            self.Integration_Result[radius] = self.Integration_Result[radius] + self.Dicom_HU_Image[
                self.Center_Row + x, self.Center_Col - y]
            if d < 0:
                d = d + 4 * x + 6
            else:
                d = d + 4 * (x - y) + 10
                y -= 1
            x += 1
    ###################################################

    def integration(self):
        for index in range(1, len(self.Integration_Result)):
            self.bresenham(index)
            self.Integration_Result[index] /= (index * 2 * 3.14)
        # calculate data by using Median
        factor = 3
        # the 1st and 2nd data = factor * md3() - md5()
        self.Median_Filter_Result[0] = np.median(self.Integration_Result[:3]) * factor - np.median(
            self.Integration_Result[:5])
        self.Median_Filter_Result[1] = np.median(self.Integration_Result[:3]) * factor - np.median(
            self.Integration_Result[:5])
        # the last and 2nd last data = factor * md3() - md5()
        self.Median_Filter_Result[-1] = np.median(self.Integration_Result[-3:]) * factor - np.median(
            self.Integration_Result[-5:])
        self.Median_Filter_Result[-2] = np.median(self.Integration_Result[-3:]) * factor - np.median(
            self.Integration_Result[-5:])
        for index in range(3, len(self.Integration_Result) - 2):
            self.Median_Filter_Result[index] = np.median(self.Integration_Result[index - 3:index + 3])
    ######################################################

    def show_image(self):
        if self.isShowImgReady:
            # set up the output file name
            try:
                image__filename = r"{0}_{1}Kv_{2}mA_{3}_{4}I{5}_{6}.jpeg".format(str(self.Dicom_Station_Name),
                                                                                 str(self.Dicom_KVP),
                                                                                 str(self.Dicom_Current),
                                                                                 str(self.Dicom_Kernel),
                                                                                 str(self.Dicom_Total_Collimation),
                                                                                 str(self.Dicom_Slice_Thickness),
                                                                                 str(self.Dicom_Instance))
                image__filename__fig = r"{0}_{1}Kv_{2}mA_{3}_{4}I{5}_{6}_fig.jpeg".format(
                    str(self.Dicom_Station_Name),
                    str(self.Dicom_KVP),
                    str(self.Dicom_Current),
                    str(self.Dicom_Kernel),
                    str(self.Dicom_Total_Collimation),
                    str(self.Dicom_Slice_Thickness),
                    str(self.Dicom_Instance))
                im = Image.fromarray(self.Image_rescale).convert("L")
            except Exception as e:
                logging.error(str(e))
                return
            # prepare to drawing the image
            draw_surface = ImageDraw.Draw(im)
            point = self.Radius
            # draw the radius circle
            bounding_box = (self.Center_Col - point, self.Center_Row - point,
                            self.Center_Col + point, self.Center_Row + point)
            draw_surface.ellipse(bounding_box)
            # save image
            try:
                im.save(self.Dicom_File_Name + image__filename, "png")
            except Exception as e:
                logging.error(str(e))
                return
            # prepare the draw the fig
            try:
                plt.plot(self.Median_Filter_Result)
                plt.savefig(self.Dicom_File_Name + image__filename__fig)
            except Exception as e:
                logging.error(str(e))
                return
            finally:
                plt.close()
        else:  # if self.isShowImgReady == False
            logging.warning(r"File is not complete initialized, skip show image.")
            return
    ######################################################

    def connect_database(self):
        if self.isShowImgReady:
            database = SQL3Handler(self.Dicom_Station_Name,
                                   self.Dicom_KVP,
                                   self.Dicom_Current,
                                   self.Dicom_Kernel,
                                   self.Dicom_Total_Collimation,
                                   self.Dicom_Slice_Thickness,
                                   self.Dicom_Instance,
                                   self.Integration_Result)
            database.insert_data()
        else:
            logging.warning(r"File is not completely initialized, skip storing in database")
            return
    # End of connect_database
    ##################################################

# End of class dicom
#######################################################################


def main():
    logging.debug(r"here is the main program")
    # if only call the script, will use ./a.dcm as input
    if len(sys.argv) == 1:
        print("Use the script as below:")
        print("python BandAssessment.py [filename]|[folder name]")
    # if with 1 parameter
    elif len(sys.argv) == 2:
        logging.debug(r"Program has 1 parameters.")
        # if a folder is given
        if os.path.isdir(sys.argv[1]):
            logging.debug(r"The parameter is a folder!")
            dicom_path = sys.argv[1]
            # the path must have a end of / or \
            if (dicom_path[-1] != "/") and (dicom_path[-1] != "\\"):
                dicom_path += "/"
                logging.warning(r"There should have a \\ at the end of path, automatically add \\")
            # list all files in the folder
            dicom_dir_list = os.listdir(dicom_path)
            for x in dicom_dir_list:
                temp = Dicom(center=0, width=100, filename=dicom_path + x)
                temp.show_image()
                temp.connect_database()
        # if a file is given
        elif os.path.isfile(sys.argv[1]):
            logging.debug(r"The parameter is a file!")
            temp = Dicom(center=0, width=100, filename=sys.argv[1])
            temp.show_image()
            temp.connect_database()
        else:
            print("Use the script as below:")
            print("python BandAssessment.py [filename]|[folder name]")


# if it is not called by a module
if __name__ == '__main__':
    main()
