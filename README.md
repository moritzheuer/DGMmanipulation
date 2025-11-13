# DGMmanipulation

Dieses Repository enthält die benötigten Algorithmen, um Digitale Geländemodelle für eine hydrologische Modellierung des Abflusses zu präparieren.
Außerdem sind die Algorithmen enthalten, um Grabenverschlüsse, Keylines und Rückhaltegräben in das DGM einzuprägen.
Die Schritte sind dabei:
  1. "DGM_preprocessing.R": Aufbereitung des DGM, Nachprägen der existierenden Gräben.
  2. "Grabenverschlüsse.R": Einprägen der Grabenverschlüsse.
  3. "Umverteilungsgräben.R": Einprägen der Umverteilungsgräben zusammen mit ihren Grabenverschlüssen und den hangabwärtigen Wällen.
  4. "Zwischengräben.R": Einprägen der Gräben zwischen den Rückegassen, zusammen mit ihren hangabwärtig liegenden Wällen.
  5. "Grabenauslässe.R": Einprägen der Wegegrabenausleitungen, zusammen mit den hangabwärtigen Wällen.
     
Enthalten sind außerdem die erstellten und verwendeten Tabellen mit den Vertexpunkten der Maßnahmen und ihren Koordinaten.
