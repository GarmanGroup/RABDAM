
$(document).ready(function() {
    /* Hides tables from view */
    $('table').hide();

    /* Underlines sublist headings when mouse is moved over them */
    $('h3').hover(function() {
        $(this).css('text-decoration', 'underline');
     }, function(){
        $(this).css('text-decoration', 'none');
      });

    /* Reveals tables when relevant sublist heading is clicked */
    $('.sublist').click(function() {
        $(this).find('table').slideToggle('slow');
    });

});
