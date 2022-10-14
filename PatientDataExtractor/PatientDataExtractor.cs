using System;
using System.Linq;
using System.Text;
using System.Collections.Generic;
using System.Reflection;
using VMS.TPS.Common.Model.API;
using VMS.TPS.Common.Model.Types;

// TODO: Replace the following version attributes by creating AssemblyInfo.cs. You can do this in the properties of the Visual Studio project.
[assembly: AssemblyVersion("1.0.0.1")]
[assembly: AssemblyFileVersion("1.0.0.1")]
[assembly: AssemblyInformationalVersion("1.0")]

// TODO: Uncomment the following line if the script requires write access.
// [assembly: ESAPIScript(IsWriteable = true)]

namespace PatientDataExtractor
{
  class Program
  {
    [STAThread]
    static void Main(string[] args)
    {
      try
      {
        using (Application app = Application.CreateApplication())
        {
          Execute(app);
        }
      }
      catch (Exception e)
      {
        Console.Error.WriteLine(e.ToString());
      }
    }
    static void Execute(Application app)
    {
            // TODO: Add your code here.
            string PatientId = "0627543";
            Patient patient = app.OpenPatientById(PatientId);
            
            IEnumerable<StructureSet> structureSets = patient.StructureSets;

            foreach (var structureSet in patient.StructureSets)
            {
                Console.Write("Structure Set ID: " + structureSet.Id + "\n");

                Console.Write("\timage origin_x : " + structureSet.Image.Origin.x + "\n");
                Console.Write("\timage origin_y : " + structureSet.Image.Origin.y + "\n");
                Console.Write("\timage origin_z : " + structureSet.Image.Origin.z + "\n");
                Console.Write("\timage user origin : " + structureSet.Image.UserOrigin + "\n");
                Console.Write("\timage res_x : " + structureSet.Image.XRes + "\n");
                Console.Write("\timage res_y : " + structureSet.Image.YRes + "\n");
                Console.Write("\timage res_z : " + structureSet.Image.ZRes + "\n");
                Console.Write("\timage size_x : " + structureSet.Image.XSize + "\n");
                Console.Write("\timage size_y : " + structureSet.Image.YSize + "\n");
                Console.Write("\timage size_z : " + structureSet.Image.ZSize + "\n");

                foreach (var image in structureSet.Image.Series.Images)
                {
                    Console.Write("\t\timage name : " + image.Id + "\n");
                }


                foreach (var structure in structureSet.Structures)
                {
                    Console.Write("\tstructure ID: " + structure.Id + "\n");
                    Console.Write("\thas segment: " + structure.HasSegment + "\n");
                    // if (!(structure.MeshGeometry is null) && name == "BODY" || name == "PTV_ALL")
                    if (structure.HasSegment && structure.Id.ToLower().Contains("body") || structure.Id == "PTV_ALL")
                    {                        
                        Console.Write("\tstructure ID: " + structure.Id + "\n");
                        Console.Write("\t\tcenter_x: " + structure.CenterPoint.x + "\n");
                        Console.Write("\t\tcenter_y: " + structure.CenterPoint.y + "\n");
                        Console.Write("\t\tcenter_z: " + structure.CenterPoint.z + "\n");
                        Console.Write("\t\tvolume: " + structure.Volume + "\n");
                        Console.Write("\t\tbounds: " + structure.MeshGeometry.Bounds + "\n");
                        Console.Write("\t\tpoints count: " + structure.MeshGeometry.Positions.Count + "\n");
                        Console.Write("\t\tindices count: " + structure.MeshGeometry.TriangleIndices.Count + "\n");
                    }
                }
            }

            foreach (var course in patient.Courses)
            {
                Console.Write("Course ID: " + course.Id + "\n");
                foreach (var setup in course.PlanSetups)
                {
                    Console.Write("\tsetup ID: " + setup.Id + "\n");
                    Console.Write("\t\tapproval status: " + setup.ApprovalStatus + "\n");                    
                    Console.Write("\t\tdose per fx: " + setup.DosePerFraction + "\n");
                    Console.Write("\t\tnumber of fx: " + setup.NumberOfFractions + "\n");
                    Console.Write("\t\ttotal dose: " + setup.TotalDose + "\n");

                    Console.Write("\t\tstructure set ID: " + setup.StructureSet.Id + "\n");

                    foreach (var structure in setup.StructureSet.Structures)
                    {                        
                        if (structure.HasSegment && structure.Id.ToLower().Contains("body") || structure.Id == "PTV_ALL")
                        {
                            Console.Write("\t\t\tstructure ID: " + structure.Id + "\n");
                            Console.Write("\t\t\t\tcenter_x: " + structure.CenterPoint.x + "\n");
                            Console.Write("\t\t\t\tcenter_y: " + structure.CenterPoint.y + "\n");
                            Console.Write("\t\t\t\tcenter_z: " + structure.CenterPoint.z + "\n");
                            Console.Write("\t\t\t\tvolume: " + structure.Volume + "\n");
                            Console.Write("\t\t\t\tbounds: " + structure.MeshGeometry.Bounds + "\n");
                            Console.Write("\t\t\t\tpoints count: " + structure.MeshGeometry.Positions.Count + "\n");
                            Console.Write("\t\t\t\tindices count: " + structure.MeshGeometry.TriangleIndices.Count + "\n");
                        }
                    }
                }
            }

                string helloworld = patient.Id;
        }
  }
}
